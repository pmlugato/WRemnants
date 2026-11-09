#!/usr/bin/env python3
"""Stage B0 plots: the kaon's Geant4 fate vs gen pT, and what it explains.

Reads the per-chunk npz files from genmatch_simtrack_fate.py and the nano-level
npz from genmatch_study_collect.py, joins them per kaon on (event, gen pt), and
writes F1-F4 plus fate_summary.md and index_fate.html into the output dir.

Classes: 0 = no end vertex in the record (the kaon left the simulated volume or
stopped by ionisation), 1 = decay in flight, 2 = hadronic inelastic with no
same-type kaon continuing, 3 = other end process, 4 = inelastic with a same-type
kaon continuing (kinked track).
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
K_EDGES = np.array([0.15, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 20.0])
R_TRK_MIN, R_TRK_GOOD = 30.0, 60.0   # cm: pixel + TIB reached; most of the tracker reached
plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.axisbelow": True, "legend.frameon": False, "figure.dpi": 130})
FIGS = []


def centers(e):
    return np.sqrt(e[:-1] * e[1:])


def xerr(e):
    c = centers(e); return np.vstack([c - e[:-1], e[1:] - c])


def binom(k, n):
    n = np.asarray(n, float); k = np.asarray(k, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(n > 0, k / n, np.nan); return p, np.where(n > 0, np.sqrt(np.clip(p * (1 - p), 0, None) / n), np.nan)


def eff(ax, pt, num, den, label, col, mk="o"):
    k = np.array([np.sum(num & den & (pt >= lo) & (pt < hi)) for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:])])
    n = np.array([np.sum(den & (pt >= lo) & (pt < hi)) for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:])])
    p, e = binom(k, n)
    ax.errorbar(centers(K_EDGES), p, yerr=e, xerr=xerr(K_EDGES), fmt=mk, color=col, ms=5, lw=1.4, label=label)
    return p


def style(ax):
    ax.set_xscale("log"); ax.set_xlim(K_EDGES[0], K_EDGES[-1]); ax.set_xlabel("gen kaon $p_T$ [GeV]"); ax.axvline(1.0, color=INK2, ls=":", lw=0.8)


def save(fig, name, cap):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight"); plt.close(fig); FIGS.append((name, cap)); print("wrote", name)


ap = argparse.ArgumentParser()
ap.add_argument("--fateDir", default="/ceph/submit/data/user/p/pmlugato/mz/studies/fate")
ap.add_argument("--nano", default="/ceph/submit/data/user/p/pmlugato/mz/studies/genmatch_study.npz")
ap.add_argument("-o", "--outdir", default=os.path.expanduser("~/public_html/mz/gen-matching-study"))
args = ap.parse_args()
OUT = args.outdir

parts = [np.load(f) for f in sorted(glob.glob(os.path.join(args.fateDir, "chunk_*.npz")))]
F = {k: np.concatenate([p[k] for p in parts]) for k in parts[0].files if not k.startswith("meta__")}
nev = int(sum(int(p["meta__nevents"]) for p in parts))
print(f"{len(parts)} chunks, {nev} events, {len(F['gen_pt'])} decays, simtrack found {F['has_simtrack'].mean():.4f}")

pt = F["gen_pt"]; acc = np.abs(F["gen_eta"]) < 2.4
cls = F["end_class"]; end_r = F["end_r"]; end_z = F["end_z"]
in_tracker = (end_r > -98) & (end_r < 110) & (np.abs(end_z) < 280)
decayed = (cls == 1) & in_tracker
interacted = (cls == 2) & in_tracker
kinked = (cls == 4) & in_tracker
other = (cls == 3) & in_tracker
# survived far enough to make a track
alive_to = lambda r: ~(in_tracker & (cls != 0) & (end_r < r))
alive30, alive60 = alive_to(R_TRK_MIN), alive_to(R_TRK_GOOD)

# ---------------------------------------------------------------- F1: fate fractions, stacked
fig, ax = plt.subplots(figsize=(9, 5.5))
cats = [("kaon reaches r > 60 cm intact", alive60, C[0]),
        ("ends between 30 and 60 cm", alive30 & ~alive60, C[3]),
        ("decays in flight before 30 cm", decayed & (end_r < R_TRK_MIN), C[1]),
        ("hadronic interaction before 30 cm (kaon gone)", interacted & (end_r < R_TRK_MIN), C[2]),
        ("inelastic before 30 cm, same-sign kaon continues (kink)", kinked & (end_r < R_TRK_MIN), C[5]),
        ("other end process before 30 cm", other & (end_r < R_TRK_MIN), C[4])]
ntot = np.array([np.sum(acc & (pt >= lo) & (pt < hi)) for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:])])
bottom = np.zeros(len(ntot))
for lab, m, col in cats:
    k = np.array([np.sum(acc & m & (pt >= lo) & (pt < hi)) for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:])])
    f = np.where(ntot > 0, k / ntot, 0)
    ax.bar(centers(K_EDGES), f, width=np.diff(K_EDGES) * 0.92, bottom=bottom, color=col, label=lab, edgecolor="white", lw=1); bottom += f
for x, n in zip(centers(K_EDGES), ntot):
    ax.text(x, 1.01, f"{n}", ha="center", fontsize=7.5, color=INK2)
style(ax); ax.set_ylim(0, 1.06); ax.set_ylabel("fraction of gen kaons (|η|<2.4)"); ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("What Geant4 does to the kaon before it can make a track", fontsize=11)
save(fig, "F1_kaon_geant_fate.png",
     "Fate of the kaon in the simulation, per gen pT, from SimTracks/SimVertices in the MC AlCaReco. Blue reaches most of the tracker; "
     "orange decayed and green interacted before the pixel+TIB region (30 cm). This is the physics loss, independent of any cut.")

# ---------------------------------------------------------------- F2: end radius
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
sel = acc & (cls == 1) & (end_r > -98)
h = axs[0].hist2d(pt[sel], np.clip(end_r[sel], 0.1, 300), bins=[np.geomspace(0.15, 20, 40), np.geomspace(0.1, 300, 60)], cmap="Blues", norm=LogNorm())
axs[0].set_yscale("log"); axs[0].set_ylabel("decay radius [cm]"); axs[0].set_title("Decay in flight: where", fontsize=11); fig.colorbar(h[3], ax=axs[0])
for r, lab in ((4.4, "pixel L1"), (30, "TIB"), (60, "TOB"), (110, "tracker edge")):
    axs[0].axhline(r, color=INK2, lw=0.8, ls="--"); axs[0].text(0.16, r * 1.1, lab, fontsize=8, color=INK2)
sel = acc & ((cls == 2) | (cls == 4)) & (end_r > -98)
h = axs[1].hist2d(pt[sel], np.clip(end_r[sel], 0.1, 300), bins=[np.geomspace(0.15, 20, 40), np.geomspace(0.1, 300, 60)], cmap="Blues", norm=LogNorm())
axs[1].set_yscale("log"); axs[1].set_ylabel("interaction radius [cm]"); axs[1].set_title("Hadronic interaction: where", fontsize=11); fig.colorbar(h[3], ax=axs[1])
for r, lab in ((4.4, "pixel L1"), (30, "TIB"), (60, "TOB"), (110, "tracker edge")):
    axs[1].axhline(r, color=INK2, lw=0.8, ls="--")
for ax in axs:
    style(ax)
save(fig, "F2_kaon_end_radius.png", "Radius at which the kaon decayed (left) or interacted (right), vs gen pT. Interactions cluster on material layers.")

# ---------------------------------------------------------------- F3: survival curves
fig, ax = plt.subplots(figsize=(9, 5.2))
eff(ax, pt, alive30, acc, "kaon intact at r = 30 cm (pixel + TIB)", C[0])
eff(ax, pt, alive60, acc, "kaon intact at r = 60 cm", C[1], "s")
eff(ax, pt, ~decayed, acc, "did not decay inside the tracker", C[2], "^")
eff(ax, pt, ~(interacted | kinked), acc, "did not interact inside the tracker", C[5], "D")
style(ax); ax.set_ylim(0.5, 1.02); ax.set_ylabel("fraction of gen kaons"); ax.legend(fontsize=9, loc="lower right")
ax.set_title("Kaon survival in the tracker vs gen pT", fontsize=11)
save(fig, "F3_kaon_survival.png", "Survival fractions. The 'intact at 30/60 cm' curves are the ceiling any tracking + selection efficiency can reach.")

# ---------------------------------------------------------------- F4: join with the nano-level study
d = np.load(args.nano)
c_ev = d["cand__event"]; c_pt = d["cand__kaon_gen_pt"]; c_trk = (d["cand__k_near_dr"] > -98) & (d["cand__k_near_dr"] < 0.1)
c_tight = (d["cand__mup_tight"] == 1) & (d["cand__mum_tight"] == 1)
c_found = d["cand__found"] == 1
key_n = {(int(e), round(float(p), 3)): i for i, (e, p) in enumerate(zip(c_ev, c_pt))}
idx = np.array([key_n.get((int(e), round(float(p), 3)), -1) for e, p in zip(F["event"], pt)])
joined = idx >= 0
print(f"joined {joined.sum()} of {len(pt)} kaons with the nano-level study")
J = joined & acc
has_trk = np.zeros(len(pt), bool); has_trk[joined] = c_trk[idx[joined]]
tight = np.zeros(len(pt), bool); tight[joined] = c_tight[idx[joined]]
found = np.zeros(len(pt), bool); found[joined] = c_found[idx[joined]]
fig, ax = plt.subplots(figsize=(9.5, 5.5))
eff(ax, pt, has_trk, J & tight, "P(kaon track | tight J/ψ)  — observed (G1)", C[1], "s")
eff(ax, pt, alive60, J & tight, "P(kaon intact at 60 cm)  — the physics ceiling", C[0])
p_cond = eff(ax, pt, has_trk, J & tight & alive60, "P(kaon track | tight J/ψ, kaon intact at 60 cm)  — tracking + producer cuts only", C[2], "^")
eff(ax, pt, has_trk, J & tight & (decayed | interacted), "P(kaon track | kaon decayed or interacted in the tracker)", C[4], "v")
style(ax); ax.set_ylim(0, 1.05); ax.set_ylabel("fraction"); ax.legend(fontsize=8.5, loc="lower right")
ax.set_title("Taking the physics out: what is left for tracking and the AlCaReco cuts", fontsize=11)
save(fig, "F4_track_given_survival.png",
     "The observed 'kaon has a track' curve (orange) against the physics ceiling (blue) and the same curve for kaons that survived intact to 60 cm (green): "
     "the green curve is the loss that tracking efficiency and the producer's cuts (DOCA, mass window) have to explain. Magenta: kaons that decayed or "
     "interacted still occasionally leave a usable track.")

# ---------------------------------------------------------------- summary
lines = [f"# Kaon fate in Geant4 (Stage B0) — {len(parts)} chunks, {nev} events, {len(pt)} exclusive decays", "",
         "| pT [GeV] | N | intact @30 cm | intact @60 cm | decayed <30 | interacted <30 | kinked <30 | observed P(track \\| tight) | P(track \\| tight, intact@60) |",
         "|---|--:|--:|--:|--:|--:|--:|--:|--:|"]
for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:]):
    b = acc & (pt >= lo) & (pt < hi); jb = b & joined & tight
    def fr(m, den):
        return f"{np.mean(m[den]):.3f}" if den.sum() else "nan"
    lines.append(f"| {lo:g}–{hi:g} | {b.sum()} | {fr(alive30, b)} | {fr(alive60, b)} | {fr(decayed & (end_r < 30), b)} | {fr(interacted & (end_r < 30), b)} | "
                 f"{fr(kinked & (end_r < 30), b)} | {fr(has_trk, jb)} | {fr(has_trk, jb & alive60)} |")
summary = "\n".join(lines)
open(os.path.join(OUT, "fate_summary.md"), "w").write(summary + "\n"); print(summary)
html = ["<html><head><meta charset='utf-8'><title>kaon fate (B0)</title><style>body{font-family:sans-serif;max-width:1400px;margin:auto} img{max-width:100%;border:1px solid #ddd} figcaption{color:#52514e;font-size:14px} pre{background:#f6f6f4;padding:12px}</style></head><body>",
        f"<h1>Stage B0 — the kaon's fate in Geant4 ({nev} MC AlCaReco events)</h1><p><a href='index.html'>back to the gen-matching study</a></p>"]
for n, c in FIGS:
    html.append(f"<figure><img src='{n}'><figcaption><b>{n}</b> — {c}</figcaption></figure>")
html.append("<pre>" + summary + "</pre></body></html>")
open(os.path.join(OUT, "index_fate.html"), "w").write("\n".join(html)); print("wrote index_fate.html")
