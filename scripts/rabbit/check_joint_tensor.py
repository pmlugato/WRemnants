#!/usr/bin/env python
"""Assert the correlation contract of a merged joint rabbit tensor.

The joint fit is only meaningful if nuisance ``A17`` is the same physical
parameter in every channel it acts on -- Z, J/psi and B+ -> J/psi K -- at the
same prefit width. Nothing in ``merge_inputs.py`` enforces that: it correlates
by name, so a mistyped map, a mismatched command-line multiplier or a flag that
was parsed but never wired produces a tensor that looks entirely healthy and
gives a wrong answer. This script is the guard. Any two or three of the channel
families may be present (Z + J/psi, Z + J/psi + B+, J/psi + B+, Z + B+).

It runs two independent families of check.

**Structural**, on the merged tensor. Every ``A{i}``/``e{i}``/``M{i}`` appears
exactly once; none acts on a single channel family only (it would correlate
nothing); no name survives with the pre-merge Z base name; the eta count implied
by the names matches what the rename map assumed.

**Provenance**, on the input histograms. The prefit widths live in the
histmaker helper, not in the tensor, so they cannot be recovered from the merged
file -- and comparing template variation sizes across channels tests nothing,
because the same shift in A changes a Z ``mll`` yield, a J/psi mass yield and a
B+ mass yield by different fractions. Instead the per-eta width of one unit of
A, e and M is computed for each input from its recorded arguments, with the
same function the helper uses (``muon_calibration.scale_param_widths``), and
compared bin by bin. Where a histmaker also recorded the widths it actually
built (``scaleVariationWidths``, recorded by the B+ histmaker), those are used,
and they must agree with the flags -- that is what catches a flag that is
recorded but never passed to the helper.

Traps this catches by construction:

* ``--scale_e`` defaults to 5.0 in the J/psi histmaker and 1.0 in mz_dilepton
* ``--calibrationStatScaling`` defaults to 2.1 in setupRabbit, Z side only
* the B+ histmaker at massfit widths against a Z/J/psi at fitted-scale widths
* the repro (scale 1/1/1) and full (5/5/1) Z/J/psi chains mixed in one merge

Exits non-zero on any failure, so a driver script can gate the fit on it.
"""

import argparse
import os
import re
import sys

import h5py
import numpy as np

from rabbit import debugdata, inputdata
from wums import ioutils

PARAMS = ["A", "e", "M"]
PARAM_RE = re.compile(r"^(A|e|M)(\d+)$")


def _s(x):
    return x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x)


class Checker:
    def __init__(self):
        self.failures = []
        self.notes = []

    def check(self, ok, message, detail=""):
        if ok:
            print(f"  PASS  {message}")
        else:
            print(f"  FAIL  {message}")
            if detail:
                print(f"        {detail}")
            self.failures.append(message)
        return ok

    def note(self, message):
        print(f"  ..    {message}")
        self.notes.append(message)


def systs_per_channel_sparse(indata, channels):
    """Which systematics touch each channel, read straight off the sparse indices.

    The obvious route -- densify and use FitDebugData.procsForNonzeroSysts --
    is not viable here. Densifying the merged tensor means materialising
    [nbinsfull, nproc, C]; on the 5D repro tensor that reached 736 GB RSS on a
    1.4 TB shared node before being killed. The sparse indices carry the same
    information for a few hundred MB.

    Layout: norm is a [nbinsfull, nproc] SparseTensor; logk is
    [nnz(norm), C] whose rows align, in canonical order, with norm's nonzero
    (bin, proc) positions. C is nsyst, or 2*nsyst when the tensor is asymmetric.
    """
    import numpy as np

    norm_idx = np.asarray(indata.norm.indices)  # [nnz, 2] -> (bin, proc)
    logk_idx = np.asarray(indata.logk.indices)  # [nnz_logk, 2] -> (row, c)
    nsyst = int(indata.nsyst)

    bins_of_row = norm_idx[:, 0]
    rows = logk_idx[:, 0]
    cols = logk_idx[:, 1]
    syst_of_entry = np.where(cols < nsyst, cols, cols - nsyst)
    bin_of_entry = bins_of_row[rows]

    out = {}
    for ch, info in indata.channel_info.items():
        if ch not in channels:
            continue
        lo, hi = int(info["start"]), int(info["stop"])
        sel = (bin_of_entry >= lo) & (bin_of_entry < hi)
        out[ch] = set(np.unique(syst_of_entry[sel]).tolist())
    return out


def read_meta_args(path, key="meta_info"):
    """Return the argparse dict recorded in a WRemnants/rabbit HDF5, or None."""
    with h5py.File(path, "r") as f:
        for k in (key, "meta", "meta_info"):
            if k not in f:
                continue
            meta = ioutils.pickle_load_h5py(f[k])
            if not isinstance(meta, dict):
                continue
            if "args" in meta:
                return meta["args"]
            # setupRabbit nests its own args one level down and the histmaker's
            # under meta_info_input.
            for sub in ("meta_info", "meta_info_input"):
                if isinstance(meta.get(sub), dict) and "args" in meta[sub]:
                    return meta[sub]["args"]
    return None


def read_setuprabbit_args(path):
    """Return setupRabbit's own args from a Z tensor (not the histmaker's)."""
    with h5py.File(path, "r") as f:
        if "meta" not in f:
            return None
        meta = ioutils.pickle_load_h5py(f["meta"])
        info = meta.get("meta_info") if isinstance(meta, dict) else None
        if isinstance(info, dict) and "args" in info:
            return info["args"]
    return None


FAMILY_LABELS = {"z": "Z", "jpsi": "J/psi", "bplus": "B+"}
# The argument that sets the number of nuisance eta groups, per histmaker.
ETA_GROUP_ARG = {"z": "etaBins", "jpsi": "etaBins", "bplus": "scaleVariationEtaBins"}
# Relative agreement required between two widths that should be identical.
WIDTH_RTOL = 1e-6


def structural_checks(c, tensor, families, z_scale_base, n_eta):
    """``families`` maps a family label to the channel names it should own."""
    print(f"\nStructural checks on {tensor}")
    indata = inputdata.FitInputData(tensor)
    is_sparse = bool(getattr(indata, "sparse", False))
    debug = None if is_sparse else debugdata.FitDebugData(indata)

    systs = [_s(s) for s in indata.systs]
    channels = list(indata.channel_info.keys())
    print(f"  channels: {channels}")
    print(f"  nsyst   : {indata.nsyst}")

    present = {}
    for fam, chans in families.items():
        got = [ch for ch in chans if ch in channels]
        c.check(
            len(got) >= 1,
            f"{fam} channel present",
            f"looked for {chans}, found {channels}",
        )
        if got:
            present[fam] = got
    c.check(
        len(present) >= 2,
        "at least two channel families to correlate",
        f"only {list(present)}",
    )

    # No pre-merge name survived the rename.
    leftover = [s for s in systs if s.startswith(z_scale_base)]
    c.check(
        not leftover,
        f"no nuisance left under the pre-merge name '{z_scale_base}*'",
        f"{len(leftover)} survived, e.g. {leftover[:5]}",
    )

    # Exactly one of each shared name.
    counts = {}
    for s in systs:
        m = PARAM_RE.match(s)
        if m:
            counts[s] = counts.get(s, 0) + 1

    expected = {f"{p}{i}" for i in range(n_eta) for p in PARAMS}
    missing = sorted(expected - set(counts))
    extra = sorted(set(counts) - expected)
    duplicated = sorted(k for k, v in counts.items() if v != 1)

    c.check(
        not missing,
        f"all {len(expected)} shared scale nuisances present",
        f"{len(missing)} missing, e.g. {missing[:5]}",
    )
    c.check(
        not extra,
        "no unexpected A/e/M names",
        f"{len(extra)} extra, e.g. {extra[:5]}",
    )
    c.check(
        not duplicated,
        "each shared nuisance appears exactly once",
        f"duplicated: {duplicated[:5]}",
    )

    # Implied eta count matches the rename map's assumption.
    indices = sorted({int(PARAM_RE.match(s).group(2)) for s in counts})
    implied = (max(indices) + 1) if indices else 0
    c.check(
        implied == n_eta,
        f"implied eta count is {n_eta}",
        f"names imply {implied} eta bins; the rename map assumed {n_eta}",
    )

    # Each shared nuisance must act on more than one family, or it is
    # correlating nothing.
    print("  checking which channel families each shared nuisance acts on ...")
    all_present = [ch for chans in present.values() for ch in chans]
    if is_sparse:
        print("  (sparse tensor: reading the index arrays, not densifying)")
        per_ch = systs_per_channel_sparse(indata, all_present)
        name_to_idx = {s: i for i, s in enumerate(systs)}

    def acts_on(name, ch):
        if is_sparse:
            return name_to_idx[name] in per_ch.get(ch, set())
        return len(debug.procsForNonzeroSysts(channels=[ch], systs=[name])) > 0

    single, partial, neither = {}, [], []
    for name in sorted(expected & set(counts), key=lambda s: (s[0], int(s[1:]))):
        fams = [
            f for f, chans in present.items() if any(acts_on(name, ch) for ch in chans)
        ]
        if len(fams) == 1:
            single.setdefault(fams[0], []).append(name)
        elif not fams:
            neither.append(name)
        elif len(fams) < len(present):
            partial.append((name, fams))

    for fam in present:
        names = single.get(fam, [])
        c.check(
            not names,
            f"no shared nuisance acts on the {fam} channels alone",
            f"{len(names)} {fam}-only, e.g. {names[:5]} -- these correlate nothing",
        )
    if partial:
        c.note(
            f"{len(partial)} shared nuisances act on only some families, e.g. "
            f"{partial[:3]}"
        )
    if neither:
        c.note(
            f"{len(neither)} shared nuisances act on no channel "
            f"(e.g. {neither[:5]}) -- trimmed as negligible, not an error"
        )

    return indata


def massfit_variances(muon_corr_data):
    """[neta, 3] massfit variances of (A, e, M) from the data correction file."""
    import ROOT

    import narf
    from wremnants.utilities import common

    path = common.calib_filepaths["data_corrfile"][muon_corr_data]
    f = ROOT.TFile.Open(path)
    out = np.stack(
        [
            narf.root_to_hist(f.Get(p), axis_names=["scale_eta"]).variances()
            for p in PARAMS
        ],
        axis=-1,
    )
    f.Close()
    return out


def implied_widths(c, label, kind, args):
    """Per fine eta bin, what one unit of (A, e, M) means in this histogram."""
    from wremnants.production import muon_calibration

    fit = bool(args.get("fitMuonScaleAndResolution"))
    if kind in ("z", "jpsi"):
        # Both histmakers pass include_covariance=not fitMuonScaleAndResolution:
        # without the flag the nuisances are eigenvectors of the massfit
        # covariance, not per-eta parameters, and A17 names nothing comparable.
        c.check(
            fit,
            f"{label}: per-eta A/e/M basis (fitMuonScaleAndResolution)",
            "without it the histmaker builds the eigen-decomposed basis, whose "
            "names do not refer to eta bins",
        )
    derived = muon_calibration.scale_param_widths(
        massfit_variances(args.get("muonCorrData")),
        scale_A=args.get("scale_A", 1.0),
        scale_e=args.get("scale_e", 1.0),
        scale_M=args.get("scale_M", 1.0),
        fit_muon_scale=fit,
        fixed_A_unc=args.get("jpsiFixedAUnc"),
        fixed_e_unc=args.get("jpsiFixedEUnc"),
        fixed_M_unc=args.get("jpsiFixedMUnc"),
    )
    recorded = args.get("scaleVariationWidths")
    if recorded is None:
        c.note(f"{label}: no recorded widths, derived from the recorded flags")
        return derived
    recorded = np.asarray(recorded, dtype=float)
    c.check(
        recorded.shape == derived.shape
        and np.allclose(recorded, derived, rtol=WIDTH_RTOL, atol=0.0),
        f"{label}: widths the helper built agree with the recorded flags",
        f"built {np.median(recorded, axis=0)} vs flags {np.median(derived, axis=0)} "
        "(eta medians) -- a flag is recorded but not reaching the helper",
    )
    return recorded


def print_width_table(widths):
    print("  one unit of (A, e, M), eta median [min, max] over the fine eta bins:")
    for label, w in widths.items():
        cells = "  ".join(
            f"{p}={np.median(w[:, i]):.3g} [{w[:, i].min():.3g}, {w[:, i].max():.3g}]"
            for i, p in enumerate(PARAMS)
        )
        print(f"    {label:<28s} {cells}")


def read_bplus_tensor_stretch(c, b_tensor, b_hist):
    """[3] factor the B+ tensor writer applied on top of the histogram widths."""
    with h5py.File(b_tensor, "r") as f:
        meta = ioutils.pickle_load_h5py(f["meta"]) if "meta" in f else {}
    base = meta.get("scale_var_scale")
    per = meta.get("scale_var_param_scale") or {}
    c.check(base is not None, f"read scale_var_scale from {b_tensor}")
    base = 1.0 if base is None else float(base)
    built_from = meta.get("input_histogram")
    if built_from is None:
        c.note(f"{b_tensor} does not record its input histogram (older tensor writer)")
    else:
        c.check(
            os.path.realpath(built_from) == os.path.realpath(b_hist),
            "the B+ tensor was built from the B+ histogram checked here",
            f"tensor input {built_from}, checked {b_hist}",
        )
    factors = np.array([base * float(per.get(p, 1.0)) for p in PARAMS])
    if np.any(factors != 1.0):
        c.note(
            "B+ tensor stretches the variations about the nominal by "
            + ", ".join(f"{p} x{f:g}" for p, f in zip(PARAMS, factors))
        )
    return factors


def provenance_checks(c, inputs, z_tensor, expect_stat_scaling, n_eta, b_tensor=None):
    """``inputs`` is a list of (label, kind, histogram path)."""
    print("\nProvenance checks on the inputs")

    all_args = {}
    for label, kind, path in inputs:
        a = read_meta_args(path)
        if c.check(a is not None, f"read recorded args from {path}"):
            all_args[label] = (kind, a)
    if len(all_args) < len(inputs):
        return

    # The nuisance eta grouping must equal the shared basis on every input.
    for label, (kind, a) in all_args.items():
        key = ETA_GROUP_ARG[kind]
        c.check(
            a.get(key) == n_eta,
            f"{label}: {key} = {n_eta}",
            f"got {a.get(key)} -- A{{i}} would group different eta ranges",
        )

    # The shared prefit widths, bin by bin. Compared against each other rather
    # than against a constant, so changing the baseline needs no edit here --
    # only that every input moved together.
    widths = {
        label: implied_widths(c, label, kind, a)
        for label, (kind, a) in all_args.items()
    }
    # The B+ tensor writer can stretch the variations about the nominal
    # (--scaleVarScale, --scaleVarParamScale); that multiplies the unit.
    if b_tensor is not None:
        b_label, b_hist = next((l, p) for l, k, p in inputs if k == "bplus")
        widths[b_label] = widths[b_label] * read_bplus_tensor_stretch(
            c, b_tensor, b_hist
        )
    print_width_table(widths)
    ref_label = next(iter(widths))
    ref = widths[ref_label]
    for label, w in widths.items():
        if label == ref_label:
            continue
        if not c.check(
            w.shape == ref.shape,
            f"{label}: same fine eta binning as {ref_label}",
            f"{w.shape} vs {ref.shape}",
        ):
            continue
        for i, p in enumerate(PARAMS):
            ratio = ref[:, i] / w[:, i]
            c.check(
                np.allclose(ratio, 1.0, rtol=WIDTH_RTOL, atol=0.0),
                f"one unit of {p} is the same in {ref_label} and {label}",
                f"{ref_label}/{label} = {np.median(ratio):.3g} (eta median), range "
                f"[{ratio.min():.3g}, {ratio.max():.3g}] -- a unit pull moves the "
                f"two channels by different shifts",
            )

    # Not a shared width, but it changes what `e` means by cosh(eta) once a leg
    # carries mass, and it is the flag that governs comparability with havyn.
    muon = {l: a for l, (k, a) in all_args.items() if k in ("z", "jpsi")}
    vals = {l: bool(a.get("nonUltraRelativisticReweight")) for l, a in muon.items()}
    if len(vals) > 1:
        c.check(
            len(set(vals.values())) == 1,
            f"nonUltraRelativisticReweight matches across the muon histmakers",
            f"{vals}",
        )
    if any(vals.values()):
        c.note(
            "nonUltraRelativisticReweight is ON: the e convention "
            "(pt vs energy, differing by cosh(eta)) now matters and this run is "
            "no longer directly comparable to havyn's -- see assumptions.md"
        )
    for l, (k, a) in all_args.items():
        if k == "bplus":
            c.note(
                f"{l}: the bachelor's energy-loss reweight uses the kaon mass "
                f"(eConvention {a.get('eConvention')}), the muon channels do not; "
                "whether e should be shared with a kaon at all is open"
            )

    if z_tensor is None:
        return
    z_hist = next((p for l, k, p in inputs if k == "z"), None)
    sr_args = read_setuprabbit_args(z_tensor)
    if c.check(sr_args is not None, f"read setupRabbit args from {z_tensor}"):
        got_inputs = sr_args.get("inputFile") or []
        c.check(
            z_hist is None
            or any(os.path.realpath(g) == os.path.realpath(z_hist) for g in got_inputs),
            "the Z tensor was built from the Z histogram checked here",
            f"tensor input {got_inputs}, checked {z_hist}",
        )
        got = sr_args.get("calibrationStatScaling")
        c.check(
            got == expect_stat_scaling,
            f"calibrationStatScaling is {expect_stat_scaling}",
            f"got {got} -- this multiplies the Z channel's scale variations "
            f"only, so the shared nuisance would be {got}x wider on the Z side",
        )
        c.check(
            bool(sr_args.get("noClosureSysts")),
            "noClosureSysts enabled on the Z tensor",
            "the closure systematics cover a correction that is no longer "
            "applied, so leaving them in double-counts the calibration",
        )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "tensor",
        nargs="?",
        default=None,
        help="merged rabbit tensor; omit with --provenanceOnly",
    )
    p.add_argument("--zHist", help="mz_dilepton histmaker HDF5")
    p.add_argument("--jpsiHist", help="dimuon_resonances_calinput histmaker HDF5")
    p.add_argument("--bHist", help="btojpsik_cvhnano histmaker HDF5")
    p.add_argument(
        "--bTensor",
        help="pre-merge B+ tensor; its recorded variation stretch multiplies the "
        "B+ widths. Without it a stretched tensor is checked as unstretched.",
    )
    p.add_argument("--zTensor", help="pre-merge Z tensor (setupRabbit)")
    p.add_argument("--zChannel", default="Z", help="post-rename Z channel")
    p.add_argument(
        "--jpsiChannels", nargs="+", default=["JPsi_prompt", "JPsi_displaced"]
    )
    p.add_argument("--bChannel", default="btojpsik_cvhnano")
    p.add_argument("--zScaleBase", default="Scale_correction_unc")
    p.add_argument("--nEta", type=int, default=48)
    p.add_argument("--calibrationStatScaling", type=float, default=1.0)
    p.add_argument(
        "--provenanceOnly",
        action="store_true",
        help="skip the structural checks on the merged tensor (cheap; this is "
        "what the fit stage re-runs as its gate)",
    )
    args = p.parse_args()

    inputs = [
        (lab, kind, path)
        for lab, kind, path in (
            ("Z", "z", args.zHist),
            ("J/psi", "jpsi", args.jpsiHist),
            ("B+", "bplus", args.bHist),
        )
        if path
    ]
    if len(inputs) < 2:
        p.error("give at least two of --zHist, --jpsiHist, --bHist")
    if args.bTensor and not args.bHist:
        p.error("--bTensor needs --bHist")
    if args.zHist and not args.zTensor:
        p.error("--zHist needs --zTensor (the Z-only setupRabbit options live there)")
    if not args.provenanceOnly and not args.tensor:
        p.error("the merged tensor is required unless --provenanceOnly")

    c = Checker()
    if not args.provenanceOnly:
        families = {}
        if args.zHist:
            families["Z"] = [args.zChannel]
        if args.jpsiHist:
            families["J/psi"] = args.jpsiChannels
        if args.bHist:
            families["B+"] = [args.bChannel]
        structural_checks(c, args.tensor, families, args.zScaleBase, args.nEta)
    provenance_checks(
        c,
        inputs,
        args.zTensor if args.zHist else None,
        args.calibrationStatScaling,
        args.nEta,
        b_tensor=args.bTensor if args.bHist else None,
    )

    print()
    if c.failures:
        print(f"FAILED -- {len(c.failures)} check(s) did not pass:")
        for f in c.failures:
            print(f"  - {f}")
        print("\nDo not fit this tensor.")
        return 1

    print(
        f"OK -- correlation contract holds, {3 * args.nEta} parameters correlated "
        f"across {', '.join(l for l, _, _ in inputs)}."
    )
    return 0


if __name__ == "__main__":
    rc = main()
    # Not sys.exit: once a ROOT file has been opened (the massfit widths), ROOT's
    # teardown in this container turns sys.exit(1) into exit status 0, and a
    # gate that cannot fail is worse than none. os._exit skips that teardown.
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(rc)
