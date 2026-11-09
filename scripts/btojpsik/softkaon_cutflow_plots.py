#!/usr/bin/env python3
"""Stage B1 plots: the cutflow from gen kaon to AlCaReco B+ candidate, per gen pT.

Reads the npz written by softkaon_cutflow.py (one row per exclusive decay) and
writes H1-H6 + cutflow_summary.md + index_cutflow.html into the study directory.
Efficiencies carry binomial errors; one y-axis per panel; fixed categorical order.
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e6e3"
K_EDGES = np.array([0.15, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5])
DOCA_CUT = 0.03
plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.axisbelow": True, "legend.frameon": False, "figure.dpi": 130})
FIGS = []


def centers(e):
    return np.sqrt(e[:-1] * e[1:])


def xerr(e):
    c = centers(e); return np.vstack([c - e[:-1], e[1:] - c])


def eff(ax, pt, num, den, label, col, mk="o", edges=K_EDGES):
    k = np.array([np.sum(num & den & (pt >= lo) & (pt < hi)) for lo, hi in zip(edges[:-1], edges[1:])], float)
    n = np.array([np.sum(den & (pt >= lo) & (pt < hi)) for lo, hi in zip(edges[:-1], edges[1:])], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(n > 0, k / n, np.nan); e = np.where(n > 0, np.sqrt(np.clip(p * (1 - p), 0, None) / n), np.nan)
    ax.errorbar(centers(edges), p, yerr=e, xerr=xerr(edges), fmt=mk, color=col, ms=5, lw=1.4, label=label)
    return p, n


def style(ax):
    ax.set_xscale("log"); ax.set_xlim(K_EDGES[0], K_EDGES[-1]); ax.set_xlabel("gen kaon $p_T$ [GeV]"); ax.axvline(1.0, color=INK2, ls=":", lw=0.8)


def save(fig, name, cap):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight"); plt.close(fig); FIGS.append((name, cap)); print("wrote", name)


ap = argparse.ArgumentParser()
ap.add_argument("-i", "--inputs", nargs="+", default=["/ceph/submit/data/user/p/pmlugato/mz/studies/softkaon/cutflow_*.npz"])
ap.add_argument("-o", "--outdir", default=os.path.expanduser("~/public_html/mz/gen-matching-study"))
args = ap.parse_args()
OUT = args.outdir
files = sorted(f for pat in args.inputs for f in glob.glob(pat))
parts = [np.load(f) for f in files]
D = {k: np.concatenate([p[k] for p in parts]) for k in parts[0].files if not k.startswith("meta__")}
nev = int(sum(int(p["meta__nevents"]) for p in parts))
pt = D["gen_pt"]; acc = np.abs(D["gen_eta"]) < 2.4
print(f"{len(parts)} files, {nev} events, {len(pt)} decays")

intact30 = ~((D["end_class"] > 0) & (D["end_r"] > -98) & (D["end_r"] < 30))
has_trk = D["trk_idx"] >= 0
good_trk = has_trk & (np.abs(D["trk_dpt"]) < 0.2)            # the kaon's own track, not a neighbour
jpsi = D["jpsi_ok"] == 1
kin = has_trk & (D["trk_pt"] > 0.1) & (np.abs(D["trk_eta"]) < 2.5)
doca_ok = jpsi & has_trk & (D["doca_max"] > -98) & (D["doca_max"] < DOCA_CUT)
mass_ok = jpsi & has_trk & (D["m_raw"] > 5.0) & (D["m_raw"] < 5.5) & (D["bpt_raw"] > 5.0)
in_std = D["in_bcand"] == 1
in_nodoca = D["in_bcand_nodoca"] == 1

# ---------------------------------------------------------------- H1: the cutflow
fig, ax = plt.subplots(figsize=(10, 6))
base = acc & jpsi
steps = [("kaon intact at 30 cm (Geant4)", intact30, C[5], "v"),
         ("a same-charge generalTrack within ΔR<0.1", has_trk, C[0], "o"),
         ("… momentum-compatible (|Δpt/pt|<0.2)  [tracking]", good_trk, C[2], "^"),
         ("… and DOCA(K, each μ) < 0.03 cm  [producer]", good_trk & doca_ok, C[1], "s"),
         ("… and m(μμK) ∈ [5.0, 5.5], B pT > 5  [producer]", good_trk & doca_ok & mass_ok, C[3], "D"),
         ("in a standard B+ candidate  [closure]", in_std, INK, "x"),
         ("in the DOCA-free B+ candidate", in_nodoca, C[4], "P")]
res = {}
for lab, m, col, mk in steps:
    res[lab] = eff(ax, pt, m, base, lab, col, mk)
style(ax); ax.set_ylim(0, 1.05); ax.set_ylabel("fraction of decays with a tight J/ψ candidate"); ax.legend(fontsize=8.5, loc="lower right")
ax.set_title("From gen kaon to AlCaReco B+ candidate, step by step (no pileup)", fontsize=11)
save(fig, "H1_cutflow_vs_pt.png",
     "Each step is cumulative, per gen kaon pT, for decays whose J/ψ was reconstructed. The gap between successive curves is that cut's loss. "
     "The black crosses must reproduce the yellow diamonds (closure of the re-implemented selection against the producer's own output).")

# ---------------------------------------------------------------- H2: step-by-step conditional efficiencies
fig, ax = plt.subplots(figsize=(10, 5.5))
eff(ax, pt, has_trk, base & intact30, "P(track | intact at 30 cm)  — tracking", C[0])
eff(ax, pt, doca_ok, base & good_trk, "P(DOCA < 0.03 | good track)  — the DOCA cut", C[1], "s")
eff(ax, pt, mass_ok, base & good_trk & doca_ok, "P(mass window & B pT | DOCA ok)", C[3], "D")
eff(ax, pt, in_std, base & good_trk & doca_ok & mass_ok, "P(standard candidate | all cuts)  — should be 1", INK, "x")
style(ax); ax.set_ylim(0, 1.05); ax.set_ylabel("conditional efficiency"); ax.legend(fontsize=8.5, loc="lower right")
ax.set_title("Each loss on its own denominator", fontsize=11)
save(fig, "H2_conditional_efficiencies.png", "The same cutflow as conditional efficiencies, so each mechanism is read directly: tracking, the DOCA cut, the kinematic window.")

# ---------------------------------------------------------------- H3: DOCA beyond the cut
fig, axs = plt.subplots(1, 2, figsize=(14, 5.2))
s = base & good_trk & (D["doca_max"] > -98)
h = axs[0].hist2d(pt[s], np.clip(D["doca_max"][s], 1e-4, 0.5), bins=[np.geomspace(0.15, 2.5, 30), np.geomspace(1e-4, 0.5, 60)], cmap="Blues", norm=LogNorm())
axs[0].set_yscale("log"); axs[0].axhline(DOCA_CUT, color=C[1], lw=1.4); axs[0].text(0.16, 0.033, "cut 0.03 cm", color=C[1], fontsize=9)
style(axs[0]); axs[0].set_ylabel("max DOCA(K track, μ tracks) [cm]"); axs[0].set_title("True kaons: DOCA including beyond the cut", fontsize=11); fig.colorbar(h[3], ax=axs[0])
for (lo, hi), col in zip([(0.15, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.5)], C):
    ss = s & (pt >= lo) & (pt < hi)
    if ss.sum():
        axs[1].hist(np.clip(D["doca_max"][ss], 0, 0.2), bins=np.linspace(0, 0.2, 81), histtype="step", lw=1.6, color=col, density=True,
                    label=f"{lo:g}–{hi:g} GeV (n={ss.sum()}, {100*np.mean(D['doca_max'][ss] < DOCA_CUT):.0f}% pass)")
axs[1].axvline(DOCA_CUT, color=INK2, ls="--"); axs[1].set_yscale("log"); axs[1].set_xlabel("max DOCA(K, μ) [cm]"); axs[1].set_ylabel("density"); axs[1].legend(fontsize=8.5)
axs[1].set_title("DOCA shape per gen pT, now with the tail the survivors could not show", fontsize=10.5)
save(fig, "H3_doca_true_kaons.png", "The producer's DOCA variable for true kaons with a good track, including the region above the cut. The legend gives the pass fraction per pT bin.")

# ---------------------------------------------------------------- H4: relaxed-cut scan
fig, ax = plt.subplots(figsize=(9, 5.2))
den = base & good_trk & (D["doca_max"] > -98)
for thr, col in zip((0.03, 0.05, 0.10, 0.20), C):
    eff(ax, pt, D["doca_max"] < thr, den, f"DOCA < {thr:g} cm", col)
sig = D["doca_max"] / np.sqrt(np.clip(D["trk_d0err"], 1e-6, None) ** 2 + 1e-6)
for thr, col, mk in ((3, C[4], "s"), (5, C[5], "D")):
    eff(ax, pt, (D["doca_max"] > -98) & (sig < thr), den, f"DOCA / σ(d0,K) < {thr}", col, mk)
style(ax); ax.set_ylim(0, 1.05); ax.set_ylabel("fraction of good kaon tracks kept"); ax.legend(fontsize=8.5, loc="lower right")
ax.set_title("What a relaxed or significance-based DOCA cut would keep", fontsize=11)
save(fig, "H4_doca_cut_scan.png", "Signal efficiency of alternative DOCA thresholds, and of a significance cut using the kaon track's σ(d0) as the DOCA error, per gen pT bin.")

# ---------------------------------------------------------------- H5: tracking loss anatomy
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
s = base & has_trk
h = axs[0].hist2d(pt[s], np.clip(D["trk_dpt"][s], -0.5, 0.5), bins=[np.geomspace(0.15, 2.5, 30), np.linspace(-0.5, 0.5, 61)], cmap="Blues", norm=LogNorm())
style(axs[0]); axs[0].set_ylabel("(pt_track − pt_gen)/pt_gen, nearest same-charge generalTrack"); axs[0].set_title("Momentum residual of the nearest generalTrack", fontsize=11); fig.colorbar(h[3], ax=axs[0])
s = base & good_trk
h = axs[1].hist2d(pt[s], D["trk_hits"][s], bins=[np.geomspace(0.15, 2.5, 30), np.arange(2.5, 30.5, 1)], cmap="Blues", norm=LogNorm())
style(axs[1]); axs[1].set_ylabel("nValidHits of the kaon generalTrack"); axs[1].set_title("Hit count of the kaon track", fontsize=11); fig.colorbar(h[3], ax=axs[1])
save(fig, "H5_tracking_anatomy.png", "Left: whether the nearest generalTrack is the kaon (core) or a neighbour (broad). Right: hit counts of the true kaon tracks.")

# ---------------------------------------------------------------- summary
lines = [f"# Stage B1 cutflow — {len(parts)} files, {nev} events, {len(pt)} exclusive decays (no PU, kaon pT < 2.5 GeV filter)", "",
         "| pT [GeV] | N (tight J/ψ) | intact@30 | track | good track | DOCA<0.03 | mass&pT | std cand | no-DOCA cand | P(DOCA \\| good trk) |", "|---|" + "--:|" * 10]
for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:]):
    b = base & (pt >= lo) & (pt < hi)
    if b.sum() == 0:
        continue
    f = lambda m: f"{np.mean(m[b]):.3f}"
    g = b & good_trk
    lines.append(f"| {lo:g}–{hi:g} | {b.sum()} | {f(intact30)} | {f(has_trk)} | {f(good_trk)} | {f(good_trk & doca_ok)} | {f(good_trk & doca_ok & mass_ok)} | {f(in_std)} | {f(in_nodoca)} | "
                 f"{np.mean(doca_ok[g]) if g.sum() else float('nan'):.3f} |")
summary = "\n".join(lines)
open(os.path.join(OUT, "cutflow_summary.md"), "w").write(summary + "\n"); print(summary)
html = ["<html><head><meta charset='utf-8'><title>soft-kaon cutflow (B1)</title><style>body{font-family:sans-serif;max-width:1400px;margin:auto} img{max-width:100%;border:1px solid #ddd} figcaption{color:#52514e;font-size:14px} pre{background:#f6f6f4;padding:12px}</style></head><body>",
        f"<h1>Stage B1 — cutflow on RECO ({nev} GEN-filtered events)</h1><p><a href='index.html'>gen-matching study</a> · <a href='index_fate.html'>B0 fate</a></p>"]
for n, c in FIGS:
    html.append(f"<figure><img src='{n}'><figcaption><b>{n}</b> — {c}</figcaption></figure>")
html.append("<pre>" + summary + "</pre></body></html>")
open(os.path.join(OUT, "index_cutflow.html"), "w").write("\n".join(html)); print("wrote index_cutflow.html")
