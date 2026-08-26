"""Loading and binning arithmetic for the CVH NanoAOD B+ -> J/psi K tensor.

The input is the `fitbins` histogram written by
`scripts/histmakers/btojpsik_cvhnano.py`:

    fit_bachelor_pt      Regular(100, 0, 10, overflow=True)
    fit_bachelor_eta     Regular(48, -2.4, 2.4)
    fit_bachelor_charge  Regular(2, -2, 2)
    jointCvh_mass        Regular(100, 5.0, 5.5, overflow=True)
    genCategory          Integer(0, 4)

with `Weight()` storage. This module reduces it to the fit binning and nothing
else -- no tensor, no fitting, no plotting -- so that the same reduction path can
later be applied to a scale-variation histogram when one exists (gap G2).

The generic helpers (`rebin_histogram`, `collapse_axes`, `assert_matching_axes`)
live in `rabbit_btojpsik_helpers.py` and are imported by the writer from there;
they are not duplicated here.
"""

import contextlib

import h5py
import hist
import numpy as np

from wremnants.production.btojpsik_cvhnano_axes import GEN_CATEGORY_LABELS
from wremnants.utilities.io_tools import base_io
from wums import logging

logger = logging.child_logger(__name__)

FITBINS = "fitbins"
PT_AXIS = "fit_bachelor_pt"
ETA_AXIS = "fit_bachelor_eta"
CHARGE_AXIS = "fit_bachelor_charge"
MASS_AXIS = "jointCvh_mass"
CATEGORY_AXIS = "genCategory"

# B+ mass, PDG. Used only to place the peak region for the seeding scale, never
# as a fit constraint.
B_MASS = 5.27934


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


@contextlib.contextmanager
def open_results(filename):
    """The histmaker output as narf's nested result dictionary, file held open.

    The dictionary holds `H5PickleProxy` objects that read on `.get()`, so the
    file has to stay open for as long as anything is being read out of it.
    Anything materialised inside the block stays valid after it; anything still
    a proxy does not.
    """
    with h5py.File(filename, "r") as h5file:
        yield base_io.load_results_h5py(h5file)


def _dataset_keys(results):
    return [
        key
        for key, value in results.items()
        if isinstance(value, dict) and "output" in value
    ]


def load_fitbins(filename, data_dataset=None, mc_dataset=None):
    """Return `(data_hist, mc_hist, resolved_names)`.

    Either dataset name may be omitted, in which case it is resolved by
    substring against the dataset keys present. Resolution is reported rather
    than done silently: reading the wrong dataset would produce a tensor that
    fits perfectly well and means something else.
    """
    with open_results(filename) as results:
        return _load_fitbins(results, filename, data_dataset, mc_dataset)


def load_named_hist(filename, dataset, name):
    """One named histogram from one dataset, materialised.

    Separate from `load_fitbins` because the variation histograms are optional
    and keyed by name: a missing one should say which name was looked for and
    what was there, not fail on a tuple unpack.
    """
    with open_results(filename) as results:
        if dataset not in results:
            raise KeyError(
                f"{filename} has no dataset {dataset!r}; found "
                f"{sorted(_dataset_keys(results))}"
            )
        output = results[dataset].get("output", {})
        if name not in output:
            raise KeyError(
                f"{filename} dataset {dataset!r} has no histogram {name!r}. "
                f"Present: {sorted(output)}. Run the histmaker with "
                "--scaleVariations."
            )
        obj = output[name]
        h = obj.get() if hasattr(obj, "get") else obj
    if h.variances() is None:
        raise RuntimeError(
            f"{name} was booked without variances, so it cannot go through "
            "reduce_fitbins. Book it with hist.storage.Weight()."
        )
    return h


def _load_fitbins(results, filename, data_dataset, mc_dataset):
    keys = _dataset_keys(results)
    if not keys:
        raise RuntimeError(f"{filename} holds no dataset with an 'output' entry")

    def resolve(requested, description):
        if requested is not None and requested in keys:
            return requested
        candidates = (
            [k for k in keys if requested.lower() in k.lower()]
            if requested is not None
            else keys
        )
        if len(candidates) != 1:
            raise RuntimeError(
                f"Could not resolve the {description} dataset from "
                f"{requested!r}: candidates {candidates or keys}"
            )
        logger.info("Resolved %s dataset to %s", description, candidates[0])
        return candidates[0]

    data_key = resolve(data_dataset, "data")
    mc_key = resolve(mc_dataset, "simulation")
    if data_key == mc_key:
        raise RuntimeError(
            f"The data and simulation datasets resolved to the same key, {data_key}"
        )

    hists = {}
    for label, key in (("data", data_key), ("simulation", mc_key)):
        output = results[key]["output"]
        if FITBINS not in output:
            raise RuntimeError(
                f"{key} has no '{FITBINS}' histogram. Re-run the histmaker with "
                f"--fitBinsHist. Available: {sorted(output.keys())[:8]}..."
            )
        # A histmaker output stores lazy proxies that read on `.get()`; a file
        # assembled from several runs stores the histograms themselves. Both are
        # legitimate inputs, and the difference is not worth making the caller
        # know about.
        stored = output[FITBINS]
        h = stored.get() if hasattr(stored, "get") else stored
        if h.variances() is None:
            raise RuntimeError(
                f"{key}'s '{FITBINS}' carries no variances, so the statistical "
                "contribution of the simulated template cannot be separated "
                "from the data's. Re-run the histmaker: fitbins must be booked "
                "with hist.storage.Weight()."
            )
        hists[label] = h

    return hists["data"], hists["simulation"], (data_key, mc_key)


# ---------------------------------------------------------------------------
# Flow handling
# ---------------------------------------------------------------------------


def _axis_index(h, name):
    return list(h.axes.name).index(name)


def fold_overflow(h, axis_name):
    """Add an axis' overflow content into its highest bin.

    `fit_bachelor_pt` stops at 10 GeV with `overflow=True`, and 24.3% of
    simulated signal sits between 5 and 10 GeV, so the tail beyond 10 is not
    empty. Dropping it would quietly remove signal from the hardest cell, which
    is the cell where the momentum scale is best measured.

    Done explicitly on values and variances rather than trusted to a rebin
    helper, and the caller is expected to check the total.
    """
    axis = h.axes[axis_name]
    if not axis.traits.overflow:
        return h

    iax = _axis_index(h, axis_name)
    view = h.view(flow=True)

    # Position of the stored bins along `iax` within the flow view.
    lo = 1 if axis.traits.underflow else 0
    last = lo + axis.size - 1

    def at(index):
        sel = [slice(None)] * view.ndim
        # Strip other axes' flow entries so the arithmetic touches only `iax`.
        for jax, other in enumerate(h.axes):
            if jax == iax:
                continue
            jlo = 1 if other.traits.underflow else 0
            sel[jax] = slice(jlo, jlo + other.size)
        sel[iax] = index
        return tuple(sel)

    out = h.copy()
    out_view = out.view(flow=False)
    over = view[at(last + 1)]
    top = view[at(last)]
    target = [slice(None)] * out_view.ndim
    target[iax] = axis.size - 1
    target = tuple(target)
    out_view["value"][target] = top["value"] + over["value"]
    out_view["variance"][target] = top["variance"] + over["variance"]
    logger.info(
        "Folded %s overflow into its highest bin: %.6g events added to %.6g",
        axis_name,
        float(np.sum(over["value"])),
        float(np.sum(top["value"])),
    )
    return out


def strip_flow(h):
    """A histogram with the same stored bins and no flow content.

    Called after `fold_overflow`, so what is discarded here is content that was
    deliberately not kept: mass outside [5.0, 5.5], which the fit window
    excludes anyway.
    """
    out = hist.Hist(*h.axes, storage=hist.storage.Weight())
    out.view(flow=False)["value"] = h.values()
    out.view(flow=False)["variance"] = h.variances()
    return out


# ---------------------------------------------------------------------------
# Binning
# ---------------------------------------------------------------------------


def _edge_index(axis, value, tol=1e-9):
    """Index of the stored edge equal to `value`, or raise.

    Raising rather than rounding is the point: a window edge that falls inside a
    stored bin would otherwise be silently moved, and the resulting tensor would
    not correspond to the window its metadata claims.
    """
    edges = np.asarray(axis.edges, dtype=float)
    matches = np.nonzero(np.isclose(edges, value, atol=tol, rtol=0.0))[0]
    if len(matches) != 1:
        raise ValueError(
            f"{value:.6g} is not an edge of axis '{axis.name}'. Its edges run "
            f"{edges[0]:.6g} to {edges[-1]:.6g} in steps of "
            f"{edges[1] - edges[0]:.6g}; the nearest edges are "
            f"{edges[np.argsort(np.abs(edges - value))[:2]]}."
        )
    return int(matches[0])


# A Gaussian's full width at half maximum, in units of sigma. The resolution
# guard below is stated on the FWHM rather than on sigma, because the pathology
# being guarded against is the *peak* fitting inside one bin, and the peak is
# about 2.4 sigma wide, not 1.
FWHM_PER_SIGMA = 2.3548200450309493

# Below this many bins per FWHM the peak is resolved, but barely. Warned about
# rather than rejected: the validated 2018 analysis ran at 1.63.
BINS_PER_FWHM_WARN = 2.0


def select_mass_window(h, low, high, nbins, sigma, axis_name=MASS_AXIS):
    """Restrict the mass axis to `[low, high]` and rebin it to `nbins` uniform bins.

    Three conditions, each reporting its numbers:

    1. both edges coincide with stored edges, so no content is split -- raises;
    2. the enclosed bin count divides by `nbins`, so no interpolation -- raises;
    3. the signal peak spans more than one bin -- raises below one bin per FWHM,
       warns below two.

    Condition 3 is the one worth having, and its threshold is the FWHM rather
    than sigma. Ten bins over the AlCaReco window [5.00, 5.50] is 50 MeV against
    a 32.5 MeV FWHM: the entire peak sits inside one bin, which removes the
    mass-scale sensitivity the channel exists to provide, and it does so while
    producing a fit that converges. That is the failure to reject. A threshold at
    sigma itself would also reject 20 MeV bins, which is what the validated 2018
    analysis used, so it would be a line drawn in the wrong place.
    """
    axis = h.axes[axis_name]
    i_lo = _edge_index(axis, low)
    i_hi = _edge_index(axis, high)
    if i_hi <= i_lo:
        raise ValueError(f"Empty mass window: {low:.6g} to {high:.6g}")

    enclosed = i_hi - i_lo
    if enclosed % nbins:
        raise ValueError(
            f"The window [{low:.6g}, {high:.6g}] encloses {enclosed} stored bins, "
            f"which does not divide by the requested {nbins}. Choose a bin count "
            f"among {sorted(n for n in range(1, enclosed + 1) if enclosed % n == 0)}."
        )
    factor = enclosed // nbins
    width = (high - low) / nbins
    fwhm = FWHM_PER_SIGMA * sigma
    bins_per_fwhm = fwhm / width
    if bins_per_fwhm <= 1.0:
        raise ValueError(
            f"Mass bin width {width * 1e3:.1f} MeV is not below the signal peak's "
            f"full width {fwhm * 1e3:.1f} MeV (sigma {sigma * 1e3:.1f} MeV), so "
            f"the whole peak would sit inside one bin and the fit would have no "
            f"mass-scale sensitivity while still converging. Narrow the window or "
            f"raise --massBins."
        )
    if bins_per_fwhm < BINS_PER_FWHM_WARN:
        logger.warning(
            "Only %.2f mass bins per signal FWHM (%.1f MeV bins, %.1f MeV FWHM). "
            "The peak is resolved but barely; the 2018 analysis ran at 1.63.",
            bins_per_fwhm,
            width * 1e3,
            fwhm * 1e3,
        )

    out = h[{axis_name: slice(i_lo, i_hi)}]
    if factor > 1:
        out = out[{axis_name: hist.rebin(factor)}]
    logger.info(
        "Mass window [%.4g, %.4g] -> %d bins of %.1f MeV "
        "(sigma %.1f MeV, %.2f bins per sigma, %.2f per FWHM)",
        low,
        high,
        nbins,
        width * 1e3,
        sigma * 1e3,
        sigma / width,
        bins_per_fwhm,
    )
    return out


def quantile_edges(h, nbins, axis_name=PT_AXIS, weights=None):
    """Quantile edges of `axis_name`, snapped to stored edges.

    Derived from whatever projection the caller passes (the simulated signal, in
    practice), because the binning has to be set by where the signal is, not by
    where the background is. Snapping means the edges are approximations of the
    true quantiles to within the stored granularity, 0.1 GeV; that is stated
    rather than hidden, and it is why the writer reports the occupancy it ends up
    with instead of the occupancy it aimed for.
    """
    axis = h.axes[axis_name]
    values = h.project(axis_name).values() if weights is None else np.asarray(weights)
    edges = np.asarray(axis.edges, dtype=float)

    total = float(values.sum())
    if total <= 0:
        raise ValueError(f"Cannot take quantiles of an empty {axis_name} projection")
    cumulative = np.concatenate([[0.0], np.cumsum(values)]) / total

    chosen = [0]
    for k in range(1, nbins):
        target = k / nbins
        # First stored edge whose cumulative fraction reaches the target, so
        # every bin holds at least its share rather than at most.
        index = int(np.searchsorted(cumulative, target, side="left"))
        index = min(max(index, chosen[-1] + 1), axis.size - (nbins - k))
        chosen.append(index)
    chosen.append(axis.size)

    if len(set(chosen)) != len(chosen):
        raise ValueError(
            f"{nbins} quantile bins cannot be placed on {axis_name} at its stored "
            f"granularity: the edge indices collided ({chosen})."
        )

    achieved = [cumulative[i] for i in chosen]
    logger.info(
        "%s quantile edges for %d bins: %s (cumulative fractions %s)",
        axis_name,
        nbins,
        [f"{edges[i]:.2f}" for i in chosen],
        [f"{a:.3f}" for a in achieved],
    )
    return [float(edges[i]) for i in chosen]


def rebin_to_edges(h, axis_name, new_edges):
    """Group an axis' stored bins into `new_edges`, which must be stored edges."""
    axis = h.axes[axis_name]
    indices = [_edge_index(axis, e) for e in new_edges]
    if indices != sorted(indices) or len(set(indices)) != len(indices):
        raise ValueError(f"{axis_name} edges are not strictly increasing: {new_edges}")

    iax = _axis_index(h, axis_name)
    values = np.add.reduceat(h.values(), indices[:-1], axis=iax)
    variances = np.add.reduceat(h.variances(), indices[:-1], axis=iax)
    # reduceat groups from each index to the next, so content beyond the last
    # requested edge is excluded by construction; take the slice explicitly so
    # a truncating edge list is visible in the totals the caller checks.
    keep = [slice(None)] * values.ndim
    keep[iax] = slice(0, len(indices) - 1)
    values = values[tuple(keep)]
    variances = variances[tuple(keep)]

    axes = list(h.axes)
    axes[iax] = hist.axis.Variable(
        new_edges, name=axis_name, underflow=False, overflow=False
    )
    out = hist.Hist(*axes, storage=hist.storage.Weight())
    out.view(flow=False)["value"] = values
    out.view(flow=False)["variance"] = variances
    return out


def select_category(h, category, axis_name=CATEGORY_AXIS):
    """One `genCategory` bin, with the axis dropped."""
    if axis_name not in h.axes.name:
        return h
    label = GEN_CATEGORY_LABELS.get(category, str(category))
    logger.info("Selecting %s category %d (%s)", axis_name, category, label)
    return h[{axis_name: hist.loc(category)}]


def sum_categories(h, categories, axis_name=CATEGORY_AXIS):
    """Several `genCategory` bins summed, with the axis dropped."""
    parts = [select_category(h, c, axis_name) for c in categories]
    out = parts[0].copy()
    for part in parts[1:]:
        out.view(flow=False)["value"] += part.values()
        out.view(flow=False)["variance"] += part.variances()
    return out


def collapse(h, axis_names):
    """Sum over the named axes, dropping them."""
    out = h
    for name in axis_names:
        if name not in out.axes.name:
            raise ValueError(f"Cannot collapse absent axis '{name}'")
        out = out[{name: hist.sum}]
        logger.info("Collapsed %s", name)
    return out


# ---------------------------------------------------------------------------
# The one reduction path
# ---------------------------------------------------------------------------


def reduce_fitbins(
    h,
    mass_low,
    mass_high,
    mass_bins,
    signal_sigma,
    pt_edges,
    eta_rebin=1,
    collapse_axes=(),
    check_totals=True,
):
    """Take a raw `fitbins` histogram to the fit binning.

    One function, applied in one order, so that a variation histogram can later
    be passed through the identical path (design decision 10) and cannot drift
    from what the nominal took.

    Order matters: the pT overflow is folded before anything else, because after
    that point the histogram is flow-free and every later step is plain
    arithmetic on stored bins.
    """
    before = float(h.values().sum())
    over = 0.0
    pt_axis = h.axes[PT_AXIS]
    if pt_axis.traits.overflow:
        iax = _axis_index(h, PT_AXIS)
        lo = 1 if pt_axis.traits.underflow else 0
        sel = [slice(None)] * h.view(flow=True).ndim
        for jax, other in enumerate(h.axes):
            jlo = 1 if other.traits.underflow else 0
            sel[jax] = slice(jlo, jlo + other.size)
        sel[iax] = lo + pt_axis.size
        over = float(np.sum(h.view(flow=True)[tuple(sel)]["value"]))

    out = fold_overflow(h, PT_AXIS)
    out = strip_flow(out)
    if check_totals:
        after_fold = float(out.values().sum())
        if not np.isclose(after_fold, before + over, rtol=1e-12, atol=1e-9):
            raise RuntimeError(
                f"Folding the {PT_AXIS} overflow did not preserve the total: "
                f"{before:.10g} + {over:.10g} != {after_fold:.10g}"
            )

    out = select_mass_window(out, mass_low, mass_high, mass_bins, signal_sigma)
    out = rebin_to_edges(out, PT_AXIS, list(pt_edges))
    if eta_rebin and eta_rebin > 1:
        if out.axes[ETA_AXIS].size % eta_rebin:
            raise ValueError(
                f"{ETA_AXIS} has {out.axes[ETA_AXIS].size} bins, which does not "
                f"divide by --etaRebin {eta_rebin}"
            )
        out = out[{ETA_AXIS: hist.rebin(eta_rebin)}]
        logger.info(
            "Rebinned %s by %d to %d bins", ETA_AXIS, eta_rebin, out.axes[ETA_AXIS].size
        )
    if collapse_axes:
        out = collapse(out, collapse_axes)
    return out


def cell_axis_names(h, exclude=(MASS_AXIS, CATEGORY_AXIS)):
    """The cell axes of a reduced histogram, in their stored order."""
    return [name for name in h.axes.name if name not in exclude]


def peak_mask(h, half_width, axis_name=MASS_AXIS, centre=B_MASS):
    """Boolean mask over `axis_name` bins whose centre lies within the peak region."""
    centres = np.asarray(h.axes[axis_name].centers, dtype=float)
    return np.abs(centres - centre) < half_width


def signal_seed_scale(data_hist, mc_hist, half_width, axis_name=MASS_AXIS):
    """Global scale putting the fitted signal normalisations near 1.

    Sideband-subtracted, and it has to be. A bare peak-region data/simulation
    ratio counts the background sitting under the peak as though it were signal:
    the data peak holds `S + B` while the template holds `S` alone, so the ratio
    comes out too big by `1 + B/S`. At the measured peak purity of roughly 20%
    that is a factor of five, and every cell then starts with a signal template
    larger than its own data -- which is exactly what a first implementation of
    this function did, leaving 24 of 36 cells with no room for any background at
    all.

    So the background under the peak is estimated from the sidebands, flat, and
    subtracted first. The estimator is crude on purpose. It is biased, and in a
    known direction: the partially reconstructed structure peaks *under* the
    signal, so a flat sideband extrapolation underestimates the background there
    and this overestimates the signal. Getting that right is the fit's job. All
    this number has to do is start the optimiser within a factor of order one,
    and the fit reports what it actually finds.
    """
    mask = peak_mask(data_hist, half_width, axis_name)
    if mask.all():
        raise ValueError(
            f"The peak region at half-width {half_width:.4g} covers the whole "
            f"window, leaving no sideband to estimate the background from. "
            f"Narrow --peakHalfWidth or widen the mass window."
        )
    iax = _axis_index(data_hist, axis_name)
    n_peak = int(mask.sum())
    n_side = int((~mask).sum())

    data_peak = float(np.sum(np.compress(mask, data_hist.values(), axis=iax)))
    data_side = float(np.sum(np.compress(~mask, data_hist.values(), axis=iax)))
    mc_peak = float(np.sum(np.compress(mask, mc_hist.values(), axis=iax)))
    if mc_peak <= 0:
        raise ValueError(
            "The simulated signal template is empty in the peak region, so no "
            "seeding scale can be derived from it."
        )

    bkg_under_peak = data_side * n_peak / n_side
    data_signal = data_peak - bkg_under_peak
    if data_signal <= 0:
        raise ValueError(
            f"Sideband subtraction leaves no signal in data: {data_peak:.6g} in "
            f"{n_peak} peak bins against an estimated {bkg_under_peak:.6g} of "
            f"background from {n_side} sideband bins. Either the window is wrong "
            f"or there is no peak to fit."
        )

    scale = data_signal / mc_peak
    logger.info(
        "Sideband-subtracted signal estimate: data peak %.6g - background %.6g "
        "= %.6g signal against %.6g simulated, giving scale %.4f (%d peak bins, "
        "%d sideband bins, implied peak purity %.3f)",
        data_peak,
        bkg_under_peak,
        data_signal,
        mc_peak,
        scale,
        n_peak,
        n_side,
        data_signal / data_peak,
    )
    return scale


def cap_seed_scale(scale, data_hist, signal_hist, headroom=0.8, quantile=5.0):
    """Cap the seed scale so cells keep room for background.

    A signal template cannot exceed the data in any cell -- the data is signal
    plus background, both non-negative -- so `scale * S_cell <= D_cell` is a
    physical bound, not a tuning knob. The bound is applied at a low percentile
    of `D/S` rather than at its strict minimum, because one nearly-empty cell
    would otherwise drive the global scale to zero; that is the only heuristic
    here, and `headroom` leaves the capped cells some background rather than
    exactly none.

    This bites in practice. The sideband estimate is biased upward, because the
    partially reconstructed structure peaks under the signal and a flat sideband
    extrapolation underestimates the background there. Measured on a 40-file
    sample the estimate came out about twice what the sparsest cells can hold:
    the data is 89% in the softest bachelor pT bin while the signal template, cut
    into pT quantiles of the signal itself, puts a third of its events in each --
    so the hard cells are where the bound is felt.
    """
    mass_index = list(signal_hist.axes.name).index(MASS_AXIS)
    signal_cell = signal_hist.values().sum(axis=mass_index).ravel()
    data_cell = data_hist.values().sum(axis=mass_index).ravel()
    populated = signal_cell > 0
    if not np.any(populated):
        return scale, None

    ratio = data_cell[populated] / signal_cell[populated]
    cap = headroom * float(np.percentile(ratio, quantile))
    if cap >= scale:
        logger.info(
            "Seed scale %.4f is within the per-cell headroom bound %.4f; kept.",
            scale,
            cap,
        )
        return scale, cap

    logger.warning(
        "Seed scale %.4f exceeds the per-cell bound %.4f and is capped to it. "
        "The sideband estimate is inconsistent with the sparsest cells, which is "
        "the expected direction: a flat extrapolation under a structure that "
        "peaks below the signal underestimates the background and so "
        "overestimates the signal. This changes only where the optimiser starts.",
        scale,
        cap,
    )
    return cap, cap


def occupancy_table(total, signal, peak_mask_1d=None, floor=100.0):
    """Per-cell occupancy summary: the criterion the binning was chosen against.

    `total` and `signal` are reduced histograms; the returned dictionary carries
    the minimum as well as the median, because the median is not what decides
    whether a cell fits.
    """
    iax = list(signal.axes.name).index(MASS_AXIS)
    signal_per_cell = signal.values().sum(axis=iax).ravel()
    total_per_cell = total.values().sum(axis=list(total.axes.name).index(MASS_AXIS))
    total_per_cell = total_per_cell.ravel()
    if peak_mask_1d is not None:
        peak_per_cell = (
            np.compress(peak_mask_1d, signal.values(), axis=iax).sum(axis=iax).ravel()
        )
    else:
        peak_per_cell = None

    def summary(array):
        if array is None:
            return None
        return {
            "sum": float(array.sum()),
            "min": float(array.min()),
            "median": float(np.median(array)),
            "max": float(array.max()),
            "cells_above_floor": int(np.sum(array >= floor)),
            "n_cells": int(array.size),
        }

    return {
        "total": summary(total_per_cell),
        "signal": summary(signal_per_cell),
        "signal_in_peak": summary(peak_per_cell),
        "floor": floor,
    }


def format_occupancy(table):
    """The occupancy table as text, for the writer's own output."""
    lines = [
        f"Per-cell occupancy (floor {table['floor']:.0f}):",
        f"{'quantity':<18}{'sum':>14}{'min':>10}{'median':>10}{'max':>12}"
        f"{'>= floor':>12}",
        "-" * 76,
    ]
    for key, label in (
        ("total", "all candidates"),
        ("signal", "signal category"),
        ("signal_in_peak", "signal in peak"),
    ):
        entry = table.get(key)
        if entry is None:
            continue
        lines.append(
            f"{label:<18}{entry['sum']:>14.6g}{entry['min']:>10.4g}"
            f"{entry['median']:>10.4g}{entry['max']:>12.6g}"
            f"{entry['cells_above_floor']:>7d}/{entry['n_cells']:<4d}"
        )
    return "\n".join(lines)
