#!/usr/bin/env python3
"""What the simulation says the non-signal mass shape actually is.

The fit currently models everything that is not the signal template with one
monotonic analytic function per cell (an exponential, or a Bernstein
polynomial). That is a modelling *choice*, and the simulation can say whether it
is the right one, because every simulated candidate carries a truth category.

The question this script answers is not "how big is the background" but "is it
monotonic". A component that peaks under the signal is not absorbable by any
monotonic function, and a fit that tries will move the signal normalisation and
the scale nuisances instead -- which is exactly the failure seen on the
observation.

Three measurements:

  1. per-category yield and shape in the fit window, weighted as the fit sees
     them;
  2. a monotonicity test per category: the best-fit exponential's residual, and
     whether the category has a local maximum inside the window;
  3. the data minus (fitted signal template) residual, against the exponential
     the fit actually uses.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_CATEGORY_LABELS,
    GEN_DATA,
    GEN_MC_CATEGORIES,
    GEN_NON_SIGNAL_CATEGORIES,
    GEN_SIGNAL,
)
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results

B_MASS = 5.27934


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-i",
        "--input",
        default="/ceph/submit/data/user/p/pmlugato/mz/histograms/btojpsik_cvhnano_fullstats.hdf5",
    )
    p.add_argument("-o", "--outdir", default=os.path.expanduser("~/public_html/mz/cvhnano-aem/bkg"))
    p.add_argument("--window", type=float, nargs=2, default=[5.05, 5.45])
    p.add_argument("--peak", type=float, nargs=2, default=[5.24, 5.32])
    return p.parse_args()


def load(filename):
    """(data mass hist, mc mass hist) projected onto (mass, genCategory)."""
    with open_results(filename) as results:
        keys = [k for k, v in results.items() if isinstance(v, dict) and "output" in v]
        data_key = next(k for k in keys if "Charmonium" in k)
        mc_key = next(k for k in keys if k != data_key)
        out = {}
        for label, key in (("data", data_key), ("mc", mc_key)):
            h = results[key]["output"]["fitbins"]
            h = h.get() if hasattr(h, "get") else h
            out[label] = h.project("jointCvh_mass", "genCategory")
        return out, data_key, mc_key


def edges_centres(h):
    ax = h.axes["jointCvh_mass"]
    return np.asarray(ax.edges), np.asarray(ax.centers)


def exp_fit(centres, values):
    """Least squares on log(values) -- the same functional form the fit uses."""
    ok = values > 0
    if ok.sum() < 3:
        return None, None
    slope, intercept = np.polyfit(centres[ok], np.log(values[ok]), 1)
    return slope, np.exp(intercept + slope * centres)


def local_max_inside(centres, values, margin=3):
    """Index of the largest interior local maximum, or None.

    Smoothed over three bins first: a single-bin statistical spike is not a
    peaking component, and calling it one would find a "peak" in every category.
    """
    kernel = np.ones(3) / 3.0
    smooth = np.convolve(values, kernel, mode="same")
    interior = np.arange(margin, len(values) - margin)
    if len(interior) == 0:
        return None
    best, best_val = None, -np.inf
    for i in interior:
        if smooth[i] >= smooth[i - 1] and smooth[i] >= smooth[i + 1] and smooth[i] > best_val:
            best, best_val = i, smooth[i]
    return best


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    hists, data_key, mc_key = load(args.input)
    lo, hi = args.window
    plo, phi = args.peak

    edges, centres = edges_centres(hists["mc"])
    inwin = (centres > lo) & (centres < hi)
    inpeak = (centres > plo) & (centres < phi)

    def cat(h, c):
        v = h[:, c].values()
        return np.asarray(v, dtype=float)

    print(f"data dataset: {data_key}")
    print(f"mc   dataset: {mc_key}")
    print(f"window [{lo}, {hi}]  peak [{plo}, {phi}]")
    print()

    mc_cats = {c: cat(hists["mc"], c) for c in GEN_MC_CATEGORIES}
    data_vals = cat(hists["data"], GEN_DATA)

    total_mc = sum(mc_cats.values())
    nonsig = sum(mc_cats[c] for c in GEN_NON_SIGNAL_CATEGORIES)

    print("=== per-category composition, simulation ===")
    print(f"{'category':14s} {'window':>12s} {'frac':>7s} {'peak':>10s} {'peakfrac':>9s} "
          f"{'expchi2/n':>10s} {'localmax':>9s}")
    rows = []
    for c in GEN_MC_CATEGORIES:
        v = mc_cats[c]
        nwin = v[inwin].sum()
        npk = v[inpeak].sum()
        slope, model = exp_fit(centres[inwin], v[inwin])
        if model is None:
            chi2n = float("nan")
        else:
            resid = v[inwin] - model
            var = np.maximum(v[inwin], 1.0)
            chi2n = float(np.sum(resid**2 / var) / max(1, inwin.sum() - 2))
        imax = local_max_inside(centres[inwin], v[inwin])
        mloc = centres[inwin][imax] if imax is not None else float("nan")
        rows.append((c, nwin, npk, chi2n, mloc, slope))
        print(f"{GEN_CATEGORY_LABELS[c]:14s} {nwin:12.1f} "
              f"{nwin/total_mc[inwin].sum():7.4f} {npk:10.1f} "
              f"{npk/max(nwin,1e-9):9.4f} {chi2n:10.2f} {mloc:9.4f}")

    print()
    print("=== the non-signal sum, which is what the analytic function stands in for ===")
    slope, model = exp_fit(centres[inwin], nonsig[inwin])
    resid = nonsig[inwin] - model
    var = np.maximum(nonsig[inwin], 1.0)
    chi2 = float(np.sum(resid**2 / var))
    ndf = int(inwin.sum() - 2)
    imax = local_max_inside(centres[inwin], nonsig[inwin])
    print(f"exponential slope      {slope:+.4f} / GeV")
    print(f"exponential chi2/ndf   {chi2:.1f} / {ndf} = {chi2/ndf:.2f}")
    print(f"interior local maximum {centres[inwin][imax]:.4f} GeV" if imax is not None
          else "interior local maximum none")
    print(f"non-signal in peak     {nonsig[inpeak].sum():.1f} of "
          f"{nonsig[inwin].sum():.1f} in window "
          f"({nonsig[inpeak].sum()/nonsig[inwin].sum():.4f})")
    print(f"signal in peak         {mc_cats[GEN_SIGNAL][inpeak].sum():.1f}")
    print(f"S/B in peak            "
          f"{mc_cats[GEN_SIGNAL][inpeak].sum()/max(nonsig[inpeak].sum(),1e-9):.3f}")

    # Which categories carry the peaking part of the non-signal sum? The
    # residual of the non-signal sum against its own exponential, decomposed.
    print()
    print("=== who carries the non-monotonic part ===")
    print("residual of the category against its OWN exponential, summed in the peak")
    for c in GEN_NON_SIGNAL_CATEGORIES:
        v = mc_cats[c][inwin]
        s, m = exp_fit(centres[inwin], v)
        if m is None:
            continue
        r = v - m
        pk = inpeak[inwin]
        print(f"{GEN_CATEGORY_LABELS[c]:14s} resid_peak {r[pk].sum():+9.1f} "
              f"of {v[pk].sum():8.1f}  ({r[pk].sum()/max(v[pk].sum(),1e-9):+.3f})")

    # ---- the data side -----------------------------------------------------
    print()
    print("=== data, and what a monotonic function can do with it ===")
    # Sideband-only exponential, extrapolated under the peak: the honest version
    # of the background the fit is asking the exponential to be.
    side = inwin & ~inpeak
    s_side, _ = exp_fit(centres[side], data_vals[side])
    ok = data_vals[side] > 0
    icept = np.polyfit(centres[side][ok], np.log(data_vals[side][ok]), 1)[1]
    model_side = np.exp(icept + s_side * centres)
    print(f"sideband exponential slope {s_side:+.4f} / GeV")
    print(f"data in window   {data_vals[inwin].sum():.0f}")
    print(f"data in peak     {data_vals[inpeak].sum():.0f}")
    print(f"bkg under peak   {model_side[inpeak].sum():.0f}")
    print(f"signal in peak   {data_vals[inpeak].sum() - model_side[inpeak].sum():.0f}")
    r_side = data_vals[side] - model_side[side]
    v_side = np.maximum(data_vals[side], 1.0)
    print(f"sideband-only exponential chi2/ndf "
          f"{np.sum(r_side**2/v_side):.1f} / {side.sum()-2} = "
          f"{np.sum(r_side**2/v_side)/(side.sum()-2):.2f}")

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), width_ratios=[1, 1])

    ax = axes[0]
    order = sorted(GEN_NON_SIGNAL_CATEGORIES, key=lambda c: -mc_cats[c][inwin].sum())
    bottom = np.zeros(inwin.sum())
    for c in order:
        v = mc_cats[c][inwin]
        ax.bar(centres[inwin], v, width=np.diff(edges)[inwin], bottom=bottom,
               label=GEN_CATEGORY_LABELS[c], linewidth=0)
        bottom += v
    ax.plot(centres[inwin], model, "k--", lw=1.8, label="exp fit to the sum")
    ax.set_yscale("log")
    ax.set_xlabel(r"$m(\mu\mu K)$ [GeV]")
    ax.set_ylabel("simulated candidates / 5 MeV")
    ax.set_xlim(lo, hi)
    ax.axvspan(plo, phi, color="0.85", zorder=0)
    ax.legend(fontsize=8, ncol=2)
    ax.set_title("non-signal categories, simulation", fontsize=10)

    ax = axes[1]
    ratio = np.divide(nonsig[inwin], model, out=np.ones_like(model), where=model > 0)
    err = np.divide(np.sqrt(np.maximum(nonsig[inwin], 0)), model,
                    out=np.zeros_like(model), where=model > 0)
    ax.errorbar(centres[inwin], ratio, yerr=err, fmt="o", ms=3, color="k")
    ax.axhline(1.0, color="r", lw=1.2)
    ax.axvspan(plo, phi, color="0.85", zorder=0)
    ax.set_xlabel(r"$m(\mu\mu K)$ [GeV]")
    ax.set_ylabel("non-signal / exponential")
    ax.set_xlim(lo, hi)
    ax.set_title(r"the monotonic model's residual, $\chi^2/\mathrm{ndf} = "
                 f"{chi2/ndf:.1f}$", fontsize=10)

    fig.tight_layout()
    out = os.path.join(args.outdir, "bkg_composition.png")
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print()
    print("wrote", out)


if __name__ == "__main__":
    main()
