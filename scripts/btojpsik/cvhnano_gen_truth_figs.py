#!/usr/bin/env python3
"""Is the signal template biased by the way the signal is defined?

The signal template is the set of candidates the producer's per-leg matcher
accepted, and the matcher requires dR < 0.02 AND |dpt|/pt < 0.1 per leg. The
second of those is a cut on momentum resolution, and momentum resolution is what
this analysis measures. So the template is selected on the quantity being
measured -- which does not by itself make it wrong, but does make it something
to measure rather than assume.

Four measurements:

  1. the gen-level B+ -> J/psi(mumu) K+ rate, found by walking the Gen table with
     no reco object involved, against the rate that reconstructs and is
     classified signal: efficiency versus kinematics, with acceptance separated
     out so the two are not multiplied together;
  2. why candidates land in the no-b-ancestor ("combinatorial") bin -- an
     unmatched leg, a failed ancestor search, or a genuinely non-b ancestor --
     and the mass shape of each, because that bin is the largest single
     component in the simulation and the fit currently replaces it with a
     monotonic function;
  3. the bachelor's relative pT residual, the quantity the matcher cut on, with
     the cut boundary drawn on it;
  4. the reconstructed bachelor spectrum of candidates called signal against the
     gen spectrum of the decays they came from.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_CATEGORY_LABELS,
    GEN_MC_CATEGORIES,
    GEN_NO_B_ANCESTOR,
    GEN_SIGNAL,
)
from wremnants.postprocessing.rabbit_cvhnano_helpers import open_results

B_MASS = 5.27934
ANC_STATE = {
    0: "a leg unmatched",
    1: "3 legs, search found nothing",
    2: "3 legs, non-b ancestor",
    3: "3 legs, b ancestor",
}
# The producer's cut, so the boundary can be drawn rather than described.
MATCH_REL_DPT = 0.1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-i", "--input",
        default="/ceph/submit/data/user/p/pmlugato/mz/histograms/"
                "btojpsik_cvhnano_gendiag_full.hdf5",
    )
    p.add_argument("-o", "--outdir",
                   default=os.path.expanduser("~/public_html/mz/cvhnano-aem/gentruth"))
    p.add_argument("--window", type=float, nargs=2, default=[5.05, 5.45])
    p.add_argument("--peak", type=float, nargs=2, default=[5.24, 5.32])
    return p.parse_args()


def load(filename):
    with open_results(filename) as results:
        keys = [k for k, v in results.items() if isinstance(v, dict) and "output" in v]
        mc = next(k for k in keys if "Charmonium" not in k)
        out = {}
        for name in ("genMatchDiag", "bachelorRelDPt", "genTruthKaon",
                     "genTruthKaonReco"):
            h = results[mc]["output"][name]
            out[name] = h.get() if hasattr(h, "get") else h
        return out, mc


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    h, mc = load(args.input)
    lo, hi = args.window
    plo, phi = args.peak
    print(f"simulation dataset: {mc}")
    print()

    # ---- 1. the gen-level denominator ------------------------------------
    truth = h["genTruthKaon"]
    reco = h["genTruthKaonReco"]
    n_all = truth.values().sum()
    n_acc = truth[{"inAcc": True}].values().sum()
    n_sig = reco[{"genCategory": GEN_SIGNAL}].values().sum()
    n_reco_any = reco.values().sum()
    print("=== 1. gen-level rate against the reconstructed signal ===")
    print(f"gen B+ -> J/psi(mumu) K+ decays                {n_all:12.1f}")
    print(f"  ... with all three legs in acceptance        {n_acc:12.1f}  "
          f"({n_acc/n_all:.4f})")
    print(f"reconstructed candidates, any category        {n_reco_any:12.1f}")
    print(f"reconstructed AND classified signal            {n_sig:12.1f}")
    print(f"efficiency wrt acceptance                                  "
          f"{n_sig/max(n_acc,1e-9):.5f}")
    print()
    print("The gap between 'in acceptance' and 'classified signal' is trigger,")
    print("reconstruction, the analysis selection AND the per-leg match. Item 2")
    print("below separates the last of these from the rest.")
    print()

    # ---- 2. why a candidate has no b ancestor ----------------------------
    diag = h["genMatchDiag"]
    mass_ax = diag.axes["jointCvh_mass"]
    centres = np.asarray(mass_ax.centers)
    inwin = (centres > lo) & (centres < hi)
    inpeak = (centres > plo) & (centres < phi)

    print("=== 2. the no-b-ancestor bin, decomposed ===")
    comb = diag[{"genCategory": GEN_NO_B_ANCESTOR}]
    tot = comb.project("jointCvh_mass").values()
    print(f"{'reason':34s} {'window':>11s} {'frac':>7s} {'peak':>10s} "
          f"{'peakfrac':>9s}")
    shapes = {}
    for state, label in ANC_STATE.items():
        v = comb[{"ancestorState": state}].project("jointCvh_mass").values()
        shapes[state] = v
        nwin, npk = v[inwin].sum(), v[inpeak].sum()
        if nwin <= 0:
            continue
        print(f"{label:34s} {nwin:11.1f} {nwin/tot[inwin].sum():7.4f} "
              f"{npk:10.1f} {npk/nwin:9.4f}")
    print(f"{'all reasons':34s} {tot[inwin].sum():11.1f} {1.0:7.4f} "
          f"{tot[inpeak].sum():10.1f} {tot[inpeak].sum()/tot[inwin].sum():9.4f}")
    print()
    # Which reason carries the local maximum at the B mass? A monotonic
    # background function cannot make one, so this decides whether the fit's
    # background model is wrong about the shape or only about the level.
    print("peak-to-flank ratio, per reason (peak mean / sideband mean).")
    print("A value above 1 is a component that PEAKS under the signal.")
    side = inwin & ~inpeak
    for state, label in ANC_STATE.items():
        v = shapes[state]
        if v[inwin].sum() <= 0:
            continue
        r = (v[inpeak].mean() / v[side].mean()) if v[side].mean() > 0 else np.inf
        print(f"  {label:34s} {r:6.3f}")
    print()

    # Same question for every category, so the combinatorial can be put in
    # context rather than judged alone.
    print("peak-to-flank ratio by truth category")
    for c in GEN_MC_CATEGORIES:
        v = diag[{"genCategory": c}].project("jointCvh_mass").values()
        if v[inwin].sum() <= 0:
            continue
        r = (v[inpeak].mean() / v[side].mean()) if v[side].mean() > 0 else np.inf
        print(f"  {GEN_CATEGORY_LABELS[c]:16s} {v[inwin].sum():10.1f}  ratio {r:7.3f}")
    print()

    # ---- 3. the quantity the matcher cut on ------------------------------
    rel = h["bachelorRelDPt"]
    rax = rel.axes["rel_dpt"]
    rc = np.asarray(rax.centers)
    print("=== 3. the bachelor relative pT residual, and the matcher's cut ===")
    print(f"the producer required |dpt|/pt < {MATCH_REL_DPT}")
    for c in (GEN_SIGNAL, GEN_NO_B_ANCESTOR):
        v = rel[{"genCategory": c}].project("rel_dpt").values()
        if v.sum() <= 0:
            continue
        tot_v = v.sum()
        mean = float(np.sum(rc * v) / tot_v)
        rms = float(np.sqrt(np.sum((rc - mean) ** 2 * v) / tot_v))
        # Population in the outermost 20% of the allowed range, i.e. pressed
        # against the cut. A flat-topped distribution that is cut off leaves a
        # large number here; a distribution comfortably inside leaves a small
        # one.
        edge = np.abs(rc) > 0.8 * MATCH_REL_DPT
        inside = np.abs(rc) < MATCH_REL_DPT
        outside = np.abs(rc) > MATCH_REL_DPT
        print(f"{GEN_CATEGORY_LABELS[c]:16s} n {tot_v:10.1f}  mean {mean:+.5f}  "
              f"rms {rms:.5f}")
        print(f"{'':16s} in the last 20% before the cut: "
              f"{v[edge].sum()/tot_v:.5f}")
        print(f"{'':16s} beyond the cut (must be 0):     "
              f"{v[outside].sum()/tot_v:.5f}")
        # How much of a Gaussian of this width the cut would remove, if the
        # width were the truth. The cut is at MATCH_REL_DPT; a wider true
        # distribution loses more.
        from math import erf, sqrt
        keep = erf(MATCH_REL_DPT / (rms * sqrt(2.0)))
        print(f"{'':16s} a Gaussian of this rms would keep {keep:.5f} "
              f"-> the cut costs {1-keep:.5f} of it")
        # And the same, split by gen pT: the resolution degrades at low pT, so
        # the cost is not uniform, and a non-uniform cost is a SHAPE bias.
        hpt = rel[{"genCategory": c}].project("gen_truth_pt", "rel_dpt")
        pt_ax = hpt.axes["gen_truth_pt"]
        print(f"{'':16s} {'gen pT':>10s} {'n':>10s} {'rms':>9s} {'cut cost':>9s}")
        for ipt in range(pt_ax.size):
            vv = hpt.values()[ipt]
            if vv.sum() < 200:
                continue
            m = float(np.sum(rc * vv) / vv.sum())
            s = float(np.sqrt(np.sum((rc - m) ** 2 * vv) / vv.sum()))
            k = erf(MATCH_REL_DPT / (s * sqrt(2.0))) if s > 0 else 1.0
            print(f"{'':16s} {pt_ax.centers[ipt]:10.2f} {vv.sum():10.1f} "
                  f"{s:9.5f} {1-k:9.5f}")
    print()

    # ---- figures ---------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))

    ax = axes[0]
    for state, label in ANC_STATE.items():
        v = shapes[state]
        if v[inwin].sum() <= 0:
            continue
        ax.step(centres[inwin], v[inwin], where="mid", lw=1.6, label=label)
    ax.step(centres[inwin], tot[inwin], where="mid", color="k", lw=1.2, ls="--",
            label="all")
    ax.axvspan(plo, phi, color="0.9", zorder=0)
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_xlabel(r"$m(\mu\mu K)$ [GeV]")
    ax.set_ylabel("simulated candidates / 5 MeV")
    ax.set_title("why a candidate has no b ancestor", fontsize=10)
    ax.legend(fontsize=8)

    ax = axes[1]
    for c, col in ((GEN_SIGNAL, "tab:blue"), (GEN_NO_B_ANCESTOR, "tab:orange")):
        v = rel[{"genCategory": c}].project("rel_dpt").values()
        if v.sum() <= 0:
            continue
        ax.step(rc, v / v.sum(), where="mid", color=col, lw=1.6,
                label=GEN_CATEGORY_LABELS[c])
    for s in (-MATCH_REL_DPT, MATCH_REL_DPT):
        ax.axvline(s, color="r", ls="--", lw=1.4)
    ax.text(MATCH_REL_DPT, ax.get_ylim()[1] * 0.6, "  matcher cut", color="r",
            fontsize=9)
    ax.set_xlabel(r"$(p_T^{\mathrm{reco}} - p_T^{\mathrm{gen}}) / p_T^{\mathrm{gen}}$,"
                  " bachelor")
    ax.set_ylabel("fraction / bin")
    ax.set_title("the matcher cuts on the resolution being measured", fontsize=10)
    ax.legend(fontsize=8)

    fig.tight_layout()
    out = os.path.join(args.outdir, "gen_truth_diag.png")
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print("wrote", out)

    # ---- 4. gen spectrum against the reconstructed one -------------------
    # The efficiency itself, per bin, because the shape overlay shows that it is
    # not flat and the direction matters: the eta groups with the largest
    # cosh(eta) lever arm on the e term are also the ones the selection depletes
    # most, so the most informative cells are the least populated.
    print()
    print("=== 4. efficiency wrt acceptance, per bin ===")
    for var, label in (("gen_truth_pt", "bachelor pT [GeV]"),
                       ("gen_truth_eta", "bachelor eta")):
        den = truth[{"inAcc": True}].project(var).values()
        num = reco[{"genCategory": GEN_SIGNAL}].project(var).values()
        cc = np.asarray(truth.axes[var].centers)
        print(f"  {label}")
        for i in range(len(cc)):
            if den[i] < 500:
                continue
            eff = num[i] / den[i]
            err = np.sqrt(max(num[i], 1.0)) / den[i]
            print(f"    {cc[i]:8.2f}  den {den[i]:9.0f}  num {num[i]:8.0f}  "
                  f"eff {eff:.5f} +- {err:.5f}")
        # And the same collapsed onto the four eta groups the fit uses, since
        # that is the binning the sensitivity is quoted in.
        if var == "gen_truth_eta":
            edges = np.asarray(truth.axes[var].edges)
            print("  the fit's four eta groups (|eta| ordering as binned)")
            n = len(cc) // 4
            for g in range(4):
                sl = slice(g * n, (g + 1) * n)
                d, u = den[sl].sum(), num[sl].sum()
                if d <= 0:
                    continue
                print(f"    group {g}  eta [{edges[g*n]:+.2f}, "
                      f"{edges[(g+1)*n]:+.2f}]  den {d:9.0f}  num {u:8.0f}  "
                      f"eff {u/d:.5f}   cosh(eta) at centre "
                      f"{np.cosh(0.5*(edges[g*n]+edges[(g+1)*n])):.3f}")

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    ax_eff = axes[2]
    for var, style in (("gen_truth_pt", None), ("gen_truth_eta", None)):
        pass
    for ax, var, label in (
        (axes[0], "gen_truth_pt", r"bachelor $p_T$ [GeV]"),
        (axes[1], "gen_truth_eta", r"bachelor $\eta$"),
    ):
        den_all = truth.project(var).values()
        den_acc = truth[{"inAcc": True}].project(var).values()
        num = reco[{"genCategory": GEN_SIGNAL}].project(var).values()
        c = np.asarray(truth.axes[var].centers)
        ax.step(c, den_all / max(den_all.sum(), 1e-9), where="mid", lw=1.6,
                label="gen, all decays")
        ax.step(c, den_acc / max(den_acc.sum(), 1e-9), where="mid", lw=1.6,
                label="gen, in acceptance")
        ax.step(c, num / max(num.sum(), 1e-9), where="mid", lw=1.6,
                label="reco, called signal")
        ax.set_xlabel(label)
        ax.set_ylabel("normalised to unit area")
        ax.legend(fontsize=8)
    axes[0].set_yscale("log")

    # Third panel: the efficiency in eta, which is what the first two imply but
    # do not state.
    den = truth[{"inAcc": True}].project("gen_truth_eta").values()
    num = reco[{"genCategory": GEN_SIGNAL}].project("gen_truth_eta").values()
    cc = np.asarray(truth.axes["gen_truth_eta"].centers)
    ok = den > 500
    ax_eff.errorbar(cc[ok], num[ok] / den[ok],
                    yerr=np.sqrt(np.maximum(num[ok], 1.0)) / den[ok],
                    fmt="o", ms=3, color="k")
    ax_eff.set_xlabel(r"bachelor $\eta$")
    ax_eff.set_ylabel("efficiency wrt acceptance")
    ax_eff.set_ylim(0.0, None)

    fig.tight_layout()
    out2 = os.path.join(args.outdir, "gen_vs_reco_spectrum.png")
    fig.savefig(out2, dpi=160, bbox_inches="tight")
    print("wrote", out2)


if __name__ == "__main__":
    main()
