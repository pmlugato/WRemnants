#!/usr/bin/env python3
"""Gen kaon pT spectrum of exclusive B+ -> J/psi(mumu) K+ decays, before any
reconstruction, from the two samples that exist (worklist section 8):

  * the MC nano Gen table (genmatch_study.npz): high statistics, but only events
    that passed the AlCaReco event gate (a reconstructed J/psi candidate);
  * the Stage B1 study re-run (softkaon/chunk_*.npz): every event the generator
    filter accepted, no reco requirement at all, kaon pT < 2.5 GeV by filter.

Both carry the production fragment's own filter (J/psi pT > 6, muon pT > 3.8,
|eta| < 2.52) and the study's kaon |eta| < 2.4 acceptance. Shapes are compared
on [0.15, 2.5] GeV; the nano is also shown out to 20 GeV, with the spectrum of
the decays that end up labelled signal in a B+ candidate, so the loss is read
directly against the unselected shape.
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e6e3"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.axisbelow": True, "legend.frameon": False, "figure.dpi": 130})
K_EDGES = np.array([0.15, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5])
LO, MID, HI = 0.15, 2.5, 20.0

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--nano", default="/ceph/submit/data/user/p/pmlugato/mz/studies/genmatch_study.npz")
ap.add_argument("--study", default="/ceph/submit/data/user/p/pmlugato/mz/studies/softkaon/chunk_*.npz")
ap.add_argument("-o", "--outdir", default=os.path.expanduser("~/public_html/mz/gen-matching-study"))
args = ap.parse_args()

# ---- nano: one row per exclusive decay found in the Gen table
N = np.load(args.nano)
n_pt, n_eta = N["K__gen_pt"], N["K__gen_eta"]
n_acc = np.abs(n_eta) < 2.4
n_sig = n_acc & (N["cand__prod_signal"] == 1)          # ends up labelled signal in a B+ candidate
n_tight = n_acc & (N["mup__muon_tight"] == 1) & (N["mum__muon_tight"] == 1)

# ---- study: one row per exclusive decay, every filtered event kept
parts = [np.load(f) for f in sorted(glob.glob(args.study))]
S = {k: np.concatenate([p[k] for p in parts]) for k in parts[0].files if not k.startswith("meta__")}
s_pt = S["gen_pt"]; s_acc = np.abs(S["gen_eta"]) < 2.4
s_tight = s_acc & (S["jpsi_ok"] == 1)
s_cand = s_acc & (S["in_bcand"] == 1)

# ---- binning: 12 log bins on [0.15, 2.5], same log step continued to ~20
step = np.log(MID / LO) / 12
edges = LO * np.exp(step * np.arange(0, 22))
edges = edges[edges < HI * 1.05]
edges_lo = edges[edges <= MID * 1.001]
ctr = np.sqrt(edges[:-1] * edges[1:]); w = np.diff(edges)


def density(pt, sel, norm_sel, edges):
    """dN/dpT per unit of the reference count (decays in norm_sel with LO<=pT<MID)."""
    n, _ = np.histogram(pt[sel], bins=edges)
    ref = np.sum(norm_sel & (pt >= LO) & (pt < MID))
    return n / ref / np.diff(edges), np.sqrt(n) / ref / np.diff(edges), n, ref


d_nano, e_nano, h_nano, ref_nano = density(n_pt, n_acc, n_acc, edges)
d_nsig, e_nsig, h_nsig, _ = density(n_pt, n_sig, n_acc, edges)
d_stud, e_stud, h_stud, ref_stud = density(s_pt, s_acc, s_acc, edges_lo)
d_scan, e_scan, h_scan, _ = density(s_pt, s_cand, s_acc, edges_lo)
nb = len(edges_lo) - 1

fig, (ax, axr) = plt.subplots(2, 1, figsize=(9.5, 7.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.15], "hspace": 0.06})
ax.errorbar(ctr, d_nano, yerr=e_nano, xerr=[ctr - edges[:-1], edges[1:] - ctr], fmt="o", ms=4.5, lw=1.4, color=C[0],
            label=f"MC nano Gen table, after the AlCaReco J/ψ gate  (N = {n_acc.sum():,})")
ax.errorbar(ctr[:nb], d_stud, yerr=e_stud, xerr=[ctr[:nb] - edges_lo[:-1], edges_lo[1:] - ctr[:nb]], fmt="s", ms=5, lw=1.4, color=C[1],
            label=f"study re-run, no reconstruction at all, gen filter K pT < 2.5  (N = {s_acc.sum():,})")
ax.step(edges, np.append(d_nsig, d_nsig[-1]), where="post", color=C[2], lw=1.8,
        label=f"nano: decays that end up labelled signal in a B+ candidate  ({n_sig.sum() / n_acc.sum():.1%} of the top curve)")
ax.fill_between(edges, np.append(d_nsig, d_nsig[-1]), step="post", color=C[2], alpha=0.12, lw=0)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(LO, HI)
ax.set_ylim(max(d_nano[d_nano > 0].min() * 0.5, 1e-4), d_nano.max() * 3)
ax.set_ylabel("dN / dp$_T$  per decay in [0.15, 2.5] GeV  [GeV$^{-1}$]")
ax.set_title("Gen kaon p$_T$ in exclusive B$^+$ → J/ψ(μμ) K$^+$ decays, kaon |η| < 2.4, before any selection", fontsize=11.5)
ax.axvline(1.0, color=INK2, ls=":", lw=0.9); ax.axvline(MID, color=INK2, ls="--", lw=0.7)
ax.text(1.0, ax.get_ylim()[1] * 0.55, " fit floor 1 GeV", color=INK2, fontsize=9, va="top")
ax.text(MID, ax.get_ylim()[1] * 0.55, " study filter 2.5 GeV", color=INK2, fontsize=9, va="top")
ax.legend(loc="lower left", fontsize=9)

r = d_stud / d_nano[:nb]; er = r * np.hypot(e_stud / d_stud, e_nano[:nb] / d_nano[:nb])
axr.errorbar(ctr[:nb], r, yerr=er, xerr=[ctr[:nb] - edges_lo[:-1], edges_lo[1:] - ctr[:nb]], fmt="s", ms=5, lw=1.4, color=C[1])
axr.axhline(1, color=INK2, lw=0.8); axr.axvline(1.0, color=INK2, ls=":", lw=0.9)
axr.set_ylim(0.4, 1.6); axr.set_ylabel("study / nano\n(shape, both to unit area\non [0.15, 2.5])", fontsize=9.5)
axr.set_xlabel("gen kaon p$_T$ [GeV]")
chi2 = np.nansum(((r - 1) / er) ** 2); ndf = np.sum(np.isfinite(er))
axr.text(0.985, 0.9, f"χ²/ndf = {chi2:.1f}/{ndf}", transform=axr.transAxes, ha="right", va="top", fontsize=9, color=INK2)
fig.savefig(os.path.join(args.outdir, "S1_gen_kaon_pt_spectrum.png"), bbox_inches="tight"); plt.close(fig)

# ---- fractions below thresholds, per sample and per selection stage
def frac_below(pt, sel, thr, lo=LO, hi=None):
    den = sel & (pt >= lo) & ((pt < hi) if hi else True); num = den & (pt < thr)
    p = num.sum() / max(den.sum(), 1); return p, np.sqrt(p * (1 - p) / max(den.sum(), 1)), int(den.sum())

rows = [("nano Gen, all in-acceptance decays (after J/ψ gate)", n_pt, n_acc),
        ("nano Gen, both muons pass the tight proxy", n_pt, n_tight),
        ("nano Gen, labelled signal in a B+ candidate", n_pt, n_sig),
        ("study re-run, all decays (no reco)", s_pt, s_acc),
        ("study re-run, tight J/ψ reconstructed", s_pt, s_tight),
        ("study re-run, in a standard B+ candidate", s_pt, s_cand)]
lines = ["# Gen kaon pT spectrum before selection (S1)", "",
         f"nano: {n_acc.sum():,} exclusive decays with kaon |eta| < 2.4 in {int(N['meta__nfiles'])} files; "
         f"study: {s_acc.sum():,} decays in {len(parts)} chunks ({int(sum(int(p['meta__nevents']) for p in parts)):,} events).", "",
         "| sample / stage | N in [0.15, 2.5] | frac. < 0.5 GeV | frac. < 1 GeV (of [0.15, 2.5]) | frac. < 1 GeV (of all, nano only) |", "|---|--:|--:|--:|--:|"]
for name, pt, sel in rows:
    p05, e05, n = frac_below(pt, sel, 0.5, hi=MID); p1, e1, _ = frac_below(pt, sel, 1.0, hi=MID)
    allp = frac_below(pt, sel, 1.0) if pt is n_pt else (np.nan, np.nan, 0)
    tail = f"{allp[0]:.3f} ± {allp[1]:.3f}" if pt is n_pt else "—"
    lines.append(f"| {name} | {n:,} | {p05:.3f} ± {e05:.3f} | {p1:.3f} ± {e1:.3f} | {tail} |")
lines += ["", "Per-bin shape ratio study / nano (both normalised to unit area on [0.15, 2.5] GeV):", "",
          "| bin [GeV] | " + " | ".join(f"{a:.2f}-{b:.2f}" for a, b in zip(edges_lo[:-1], edges_lo[1:])) + " |",
          "|---|" + "--:|" * nb,
          "| ratio | " + " | ".join(f"{x:.2f}±{y:.2f}" for x, y in zip(r, er)) + " |",
          "", f"chi2/ndf of ratio = 1: {chi2:.1f}/{ndf}"]
open(os.path.join(args.outdir, "gen_spectrum_summary.md"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("nano min gen pt:", n_pt[n_acc].min(), " study min gen pt:", s_pt[s_acc].min())
