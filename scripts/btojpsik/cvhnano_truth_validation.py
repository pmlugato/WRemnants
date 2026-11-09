"""Validate the regenerated CVH NanoAOD against what it was asked to fix.

Every number this prints answers one closed gap, and the point of collecting
them in one script is that a gap closed on paper that does not move its
disagreement is a finding, not a success. Reported:

  * the truth composition, by category, over the full range and in the peak
  * how much of each category lands in the peak window -- the test that
    separates a signal template from a partially reconstructed decay
  * the kaon momentum response (gap G2, which previously had no input at all)
  * data/simulation agreement on the impact parameter, before the fix
    (origin-referenced `dxy`) and after (PV-referenced `d0`), against the
    baseline the downstream workaround reached
  * data/simulation agreement on the pileup observables under true-interaction
    reweighting, against the baseline the `nPV` stopgap reached

    python scripts/btojpsik/cvhnano_truth_validation.py <histmaker.hdf5>
"""

import argparse
import json

import numpy as np

from wremnants.postprocessing import rabbit_cvhnano_helpers as H
from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_CATEGORY_LABELS,
    GEN_DATA,
    GEN_MC_CATEGORIES,
    GEN_SIGNAL,
    GEN_TEMPLATE_CATEGORIES,
)

DATA = "Charmonium_2016PostVFP_cvhNano"
MC = "InclusiveBToJpsiX_2016PostVFP"
PEAK = (5.24, 5.32)

# Baselines this run has to reproduce or beat, from Phase 2 of
# add-jpsix-histmaker. Quoted here so a regression is visible in the output
# rather than in someone's memory.
BASELINE = {
    "nPV": "201.4 -> 1.34 under the nPV stopgap",
    "impact parameter": "2340 -> 3.92 under the downstream correction",
    "nValidHits": "112.9 -> 110.3, i.e. unmoved",
}


def get(o, name):
    h = o[name]
    return h.get() if hasattr(h, "get") else h


def chi2_per_bin(d, dv, m, mv, min_entries=20.0):
    """Shape chi-square per bin between two histograms, each normalised to unit area.

    Normalised, because the simulation's cross section is a placeholder, so only
    the shapes are comparable. Bins where either sample has fewer than
    `min_entries` are dropped: a chi-square contribution from three events is
    noise that dominates the sum without meaning anything.
    """
    keep = (d > min_entries) & (m > min_entries)
    if keep.sum() < 3:
        return float("nan"), 0
    dn, mn = d / d.sum(), m / m.sum()
    dnv, mnv = dv / d.sum() ** 2, mv / m.sum() ** 2
    num = (dn[keep] - mn[keep]) ** 2
    den = dnv[keep] + mnv[keep]
    return float(np.sum(num / den) / keep.sum()), int(keep.sum())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("infile")
    p.add_argument("-o", "--output", default=None, help="JSON summary")
    p.add_argument("--arm", default="jointCvh")
    args = p.parse_args()

    out = {}
    with H.open_results(args.infile) as r:
        mo = r[MC]["output"]
        do = r[DATA]["output"] if DATA in r else None

        # --- composition ---------------------------------------------------
        hm = get(mo, f"nominal_BuJpsiK_{args.arm}Mass")
        axis = hm.axes[0]
        peak = (axis.centers > PEAK[0]) & (axis.centers < PEAK[1])
        v = hm.values()
        comp = {}
        peak_total = v[peak][:, GEN_MC_CATEGORIES].sum()
        for c in GEN_MC_CATEGORIES:
            n_all = float(v[:, c].sum())
            n_peak = float(v[peak, c].sum())
            comp[GEN_CATEGORY_LABELS[c]] = dict(
                total=n_all,
                in_peak=n_peak,
                frac_of_category_in_peak=n_peak / n_all if n_all > 0 else 0.0,
                pct_of_peak=100.0 * n_peak / peak_total if peak_total > 0 else 0.0,
            )
        out["composition"] = comp
        real_b_bkg = sum(
            comp[GEN_CATEGORY_LABELS[c]]["in_peak"]
            for c in GEN_TEMPLATE_CATEGORIES
            if c != GEN_SIGNAL
        )
        sig_peak = comp[GEN_CATEGORY_LABELS[GEN_SIGNAL]]["in_peak"]
        out["real_b_background_over_signal"] = (
            real_b_bkg / sig_peak if sig_peak > 0 else float("nan")
        )

        print(f"Truth composition, {args.arm} mass, peak window {PEAK}")
        print(f"{'category':14s}{'total':>12s}{'in peak':>11s}{'in-peak frac':>14s}{'% of peak':>11s}")
        for name, c in comp.items():
            print(f"{name:14s}{c['total']:12.1f}{c['in_peak']:11.1f}"
                  f"{c['frac_of_category_in_peak']:14.3f}{c['pct_of_peak']:11.2f}")
        print(f"\nreal-b background / signal in the peak: "
              f"{out['real_b_background_over_signal']:.4f}")

        # --- kaon momentum response (gap G2) -------------------------------
        if "nominal_cand_bachelorPtResponse" in mo:
            hr = get(mo, "nominal_cand_bachelorPtResponse")
            rv = hr.values()[:, GEN_SIGNAL]
            ce = hr.axes[0].centers
            if rv.sum() > 0:
                mean = float((ce * rv).sum() / rv.sum())
                rms = float(np.sqrt((rv * (ce - mean) ** 2).sum() / rv.sum()))
                out["kaon_response"] = dict(mean=mean, rms=rms, entries=float(rv.sum()))
                print(f"\nkaon pT response, signal category: mean {mean:.5f}, "
                      f"rms {rms:.4f}, {rv.sum():.0f} candidates")

        if do is None:
            print("\nno data in this file; skipping the data/simulation tests")
            if args.output:
                json.dump(out, open(args.output, "w"), indent=1)
            return

        # --- data/simulation shape tests -----------------------------------
        tests = [
            ("bachelor dxy (origin-referenced, the defect)", "nominal_cand_bachelorDxy"),
            ("bachelor d0 (PV-referenced, the fix)", "nominal_cand_bachelorD0"),
            ("bachelor d0 significance", "nominal_cand_bachelorD0Sig"),
            ("bachelor nValidHits", "nominal_cand_bachelorNValidHits"),
            ("nPV", "nominal_nPV"),
            (f"{args.arm} mass", f"nominal_BuJpsiK_{args.arm}Mass"),
            ("kaon pT", "nominal_BuJpsiK_kaonPt"),
            ("unconstrained dimuon mass", "nominal_BuJpsiK_kvfCvhDimuonMass"),
        ]
        print(f"\n{'observable':46s}{'chi2/bin':>10s}{'bins':>6s}")
        out["shape_tests"] = {}
        for label, hname in tests:
            if hname not in do or hname not in mo:
                continue
            hd, hmm = get(do, hname), get(mo, hname)
            d = hd.values()[:, GEN_DATA]
            dv = hd.variances()[:, GEN_DATA]
            m = hmm.values()[:, GEN_MC_CATEGORIES].sum(axis=1)
            mv = hmm.variances()[:, GEN_MC_CATEGORIES].sum(axis=1)
            c2, nb = chi2_per_bin(d, dv, m, mv)
            out["shape_tests"][label] = dict(chi2_per_bin=c2, bins=nb)
            print(f"{label:46s}{c2:10.2f}{nb:6d}")

        # Fraction within 0.05 cm of the vertex, the number the workaround
        # reached 0.817 (data) / 0.844 (simulation) on.
        for label, hname in [("dxy", "nominal_cand_bachelorDxy"),
                             ("d0", "nominal_cand_bachelorD0")]:
            if hname not in do:
                continue
            hd, hmm = get(do, hname), get(mo, hname)
            ce = hd.axes[0].centers
            near = np.abs(ce) < 0.05
            fd = float(hd.values()[near, GEN_DATA].sum() / hd.values()[:, GEN_DATA].sum())
            mvals = hmm.values()[:, GEN_MC_CATEGORIES].sum(axis=1)
            fm = float(mvals[near].sum() / mvals.sum())
            out.setdefault("within_0p05cm", {})[label] = dict(data=fd, mc=fm)
            print(f"|{label}| < 0.05 cm: data {fd:.3f}, simulation {fm:.3f}")

        print("\nBaselines to beat:")
        for k, v_ in BASELINE.items():
            print(f"  {k}: {v_}")

    if args.output:
        json.dump(out, open(args.output, "w"), indent=1)
        print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
