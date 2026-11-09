#!/usr/bin/env python3
"""Gen-matching study, plotting step (worklist 3.1-3.4).

Reads the .npz written by genmatch_study_collect.py and produces the figures
plus a summary table and an index.html into the output directory.

Conventions: points with binomial errors for every efficiency, ratio panels on
every A/B shape comparison, one y-axis per panel, categorical colours in a fixed
order (blue, orange, aqua, yellow, magenta, violet), sequential single-hue
colormap for 2D densities. Counts are unweighted (see the collector).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# Validated categorical order (dataviz palette.md, light surface).
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]
INK = "#0b0b0b"; INK2 = "#52514e"; GRID = "#e6e6e3"
PROD_DR, PROD_DPT = 0.02, 0.1
SENT = -98.0

# Gen kaon pT bins (GeV). The fit selection starts at 1.0; the soft region is
# resolved finely because that is where the question lives.
K_EDGES = np.array([0.15, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 20.0])
MU_EDGES = np.array([2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 25.0])
COARSE = [(0.15, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 5.0), (5.0, 20.0)]

plt.rcParams.update({
    "font.size": 11, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK,
    "ytick.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "legend.frameon": False, "figure.dpi": 130,
})

FIGS = []  # (filename, caption)


def save(fig, name, caption):
    path = os.path.join(OUT, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    FIGS.append((name, caption))
    print("wrote", name)


def binom(k, n):
    n = np.asarray(n, float); k = np.asarray(k, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(n > 0, k / n, np.nan)
        e = np.where(n > 0, np.sqrt(np.clip(p * (1 - p), 0, None) / n), np.nan)
    return p, e


def centers(edges):
    return np.sqrt(edges[:-1] * edges[1:])  # geometric, for log-x


def xerr(edges):
    c = centers(edges)
    return np.vstack([c - edges[:-1], edges[1:] - c])


def eff_vs_pt(ax, pt, num, den, edges, label, color, marker="o"):
    """Draw efficiency num/den in bins of pt with binomial errors."""
    k = np.array([np.sum(num & den & (pt >= lo) & (pt < hi)) for lo, hi in zip(edges[:-1], edges[1:])])
    n = np.array([np.sum(den & (pt >= lo) & (pt < hi)) for lo, hi in zip(edges[:-1], edges[1:])])
    p, e = binom(k, n)
    ax.errorbar(centers(edges), p, yerr=e, xerr=xerr(edges), fmt=marker, color=color, ms=5,
                lw=1.4, capsize=0, label=label)
    return p, e, n


def style_pt_axis(ax, edges, xlabel="gen kaon $p_T$ [GeV]"):
    ax.set_xscale("log"); ax.set_xlim(edges[0], edges[-1]); ax.set_xlabel(xlabel)
    ax.axvline(1.0, color=INK2, lw=0.8, ls=":", zorder=0)


def quantile_curves(ax, x, y, edges, qs=(2.5, 16, 50, 84, 97.5), color=INK):
    cx = centers(edges)
    for q, ls in zip(qs, [":", "--", "-", "--", ":"]):
        vals = [np.percentile(y[(x >= lo) & (x < hi)], q) if np.sum((x >= lo) & (x < hi)) > 20 else np.nan
                for lo, hi in zip(edges[:-1], edges[1:])]
        ax.plot(cx, vals, color=color, ls=ls, lw=1.2, label=f"p{q:g}" if q in (50, 84, 97.5) else None)


# ---------------------------------------------------------------- loading
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", default="/ceph/submit/data/user/p/pmlugato/mz/studies/genmatch_study.npz")
ap.add_argument("-o", "--outdir", default=os.path.expanduser("~/public_html/mz/gen-matching-study"))
args = ap.parse_args()
OUT = args.outdir
os.makedirs(OUT, exist_ok=True)
d = np.load(args.input)
nfiles = int(d["meta__nfiles"]); nev = int(d["meta__nevents"])

K = {k.split("__", 1)[1]: d[k] for k in d.files if k.startswith("K__")}
MU = {}
for k in d.files:
    if k.startswith("mup__") or k.startswith("mum__"):
        MU.setdefault(k.split("__", 1)[1], []).append(d[k])
MU = {k: np.concatenate(v) for k, v in MU.items()}
CD = {k.split("__", 1)[1]: d[k] for k in d.files if k.startswith("cand__")}

nK = len(K["gen_pt"])
print(f"{nfiles} files, {nev} events, {nK} exclusive B+ -> J/psi(mumu) K decays")

# Acceptance: the kaon must be where a track can exist at all.
acc = np.abs(K["gen_eta"]) < 2.4
kpt = K["gen_pt"]
has_trk = K["near_dr"] > SENT
in_loose = K["n_in_loose"] > 0            # a same-charge track within dR < 0.1
near_dr_ok = has_trk & (K["near_dr"] < PROD_DR)
near_dpt_ok = has_trk & (np.abs(K["near_dpt"]) < PROD_DPT)
anypass = K["any_pass"] == 1
# best |dpt| among tracks inside the dR cone, and best dR among pt-compatible tracks
bdpt = K["best_dpt_within_dr"]; has_bdpt = bdpt > SENT
bdr = K["best_dr_passing_dpt"]; has_bdr = bdr > SENT

# ================================================================ K01: dpt vs pt, 2D
sel = acc & near_dr_ok
fig, ax = plt.subplots(figsize=(7.5, 5))
h = ax.hist2d(kpt[sel], K["near_dpt"][sel], bins=[np.geomspace(0.15, 20, 60), np.linspace(-0.4, 0.4, 81)],
              cmap="Blues", norm=LogNorm())
quantile_curves(ax, kpt[sel], K["near_dpt"][sel], K_EDGES)
for y in (-PROD_DPT, PROD_DPT):
    ax.axhline(y, color=C[1], lw=1.4)
ax.text(0.16, PROD_DPT + 0.01, "producer cut |Δpt/pt| < 0.1", color=C[1], fontsize=9)
style_pt_axis(ax, np.array([0.15, 20]))
ax.set_ylabel("(pt_reco − pt_gen) / pt_gen, nearest same-charge track (dR < 0.02)")
ax.set_title("Kaon momentum residual vs gen pT  (curves: p2.5/16/50/84/97.5)", fontsize=11)
fig.colorbar(h[3], ax=ax, label="decays / bin")
ax.legend(loc="lower right", fontsize=9)
save(fig, "K01_kaon_dpt_vs_pt_2d.png",
     "Kaon Δpt/pt of the nearest same-charge track within dR<0.02, against gen kaon pT. Orange: the producer's ±0.1 window. "
     "Quantile curves show where the residual tails sit relative to the cut.")

# ================================================================ K02: dR vs pt, 2D (two panels)
fig, axs = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
for ax, s, title in ((axs[0], acc & has_trk, "nearest same-charge track, any pt"),
                     (axs[1], acc & has_bdr & in_loose, "nearest track with |Δpt/pt| < 0.1 (a track within 0.1 exists)")):
    y = K["near_dr"][s] if ax is axs[0] else bdr[s]
    h = ax.hist2d(kpt[s], np.clip(y, 1e-4, 2), bins=[np.geomspace(0.15, 20, 60), np.geomspace(1e-4, 2, 70)],
                  cmap="Blues", norm=LogNorm())
    ax.set_yscale("log")
    quantile_curves(ax, kpt[s], np.clip(y, 1e-4, 2), K_EDGES, qs=(50, 84, 97.5))
    ax.axhline(PROD_DR, color=C[1], lw=1.4); ax.text(0.16, PROD_DR * 1.15, "dR < 0.02", color=C[1], fontsize=9)
    style_pt_axis(ax, np.array([0.15, 20])); ax.set_title(title, fontsize=11)
axs[0].set_ylabel("ΔR(gen kaon, reco track)"); axs[0].legend(loc="upper right", fontsize=9)
save(fig, "K02_kaon_dR_vs_pt_2d.png",
     "ΔR between the gen kaon and reco tracks vs gen pT. Left: nearest same-charge track regardless of momentum (includes events where the kaon "
     "was not reconstructed and the nearest track is unrelated). Right: nearest momentum-compatible track, given some track exists within 0.1. Entries far above 0.1 are unrelated tracks: the kaon itself is not in the Track table (it entered no AlCaReco candidate), so no threshold choice can recover them. Orange: producer cut.")

# ================================================================ K03: breakdown, stacked, denominator = acceptance
cats = [("pass (dR<0.02 & |Δpt/pt|<0.1, any track)", anypass),
        ("nearest fails dR only", ~anypass & has_trk & (K["near_dr"] >= PROD_DR) & in_loose & near_dpt_ok),
        ("nearest fails |Δpt/pt| only", ~anypass & near_dr_ok & ~near_dpt_ok),
        ("nearest fails both (track within 0.1)", ~anypass & in_loose & ~near_dr_ok & ~near_dpt_ok),
        ("no same-charge track within ΔR<0.1", ~in_loose)]
fig, ax = plt.subplots(figsize=(8.5, 5.2))
bottom = np.zeros(len(K_EDGES) - 1)
ntot = np.array([np.sum(acc & (kpt >= lo) & (kpt < hi)) for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:])])
for (lab, m), col in zip(cats, C):
    k = np.array([np.sum(acc & m & (kpt >= lo) & (kpt < hi)) for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:])])
    f = np.where(ntot > 0, k / ntot, 0)
    ax.bar(centers(K_EDGES), f, width=np.diff(K_EDGES) * 0.92, bottom=bottom, color=col, label=lab, edgecolor="white", lw=1)
    bottom += f
for x, n in zip(centers(K_EDGES), ntot):
    ax.text(x, 1.01, f"{n}", ha="center", fontsize=7.5, color=INK2)
style_pt_axis(ax, K_EDGES); ax.set_ylim(0, 1.06); ax.set_ylabel("fraction of gen kaons (|η|<2.4)")
ax.set_title("Where the kaon match is lost, per gen pT bin (numbers: decays per bin)", fontsize=11)
ax.legend(loc="lower left", fontsize=8.5)
save(fig, "K03_kaon_match_breakdown.png",
     "Fate of every gen kaon in acceptance: matched by the producer's criteria, or lost to the angular cut, to the momentum cut, to both, "
     "or because no same-charge track exists within ΔR<0.1 (not reconstructed / not in the AlCaReco track collection).")

# ================================================================ K04: efficiencies given a track exists
fig, ax = plt.subplots(figsize=(8.5, 5.2))
den = acc & in_loose
eff_vs_pt(ax, kpt, anypass, den, K_EDGES, "P(pass both | track within 0.1)", C[0])
eff_vs_pt(ax, kpt, near_dr_ok, den, K_EDGES, "P(nearest track has dR<0.02, any pt | track)  — the angular cut alone", C[1], "s")
eff_vs_pt(ax, kpt, has_bdpt & (np.abs(bdpt) < PROD_DPT), den, K_EDGES, "P(some track in dR<0.02 has |Δpt/pt|<0.1 | track)", C[2], "^")
eff_vs_pt(ax, kpt, has_bdpt & (np.abs(bdpt) < PROD_DPT), acc & has_bdpt, K_EDGES, "P(|Δpt/pt|<0.1 | dR<0.02)  — the pt cut's own bite", C[5], "D")
style_pt_axis(ax, K_EDGES); ax.set_ylim(0, 1.05); ax.set_ylabel("efficiency"); ax.legend(loc="lower right", fontsize=8.5)
ax.set_title("Kaon matching efficiency conditional on a reconstructed track", fontsize=11)
save(fig, "K04_kaon_eff_given_track.png",
     "The matching question proper: given a same-charge track within ΔR<0.1, how often does the producer's pair of cuts accept it, and which cut bites. "
     "Violet: the momentum cut's efficiency for tracks that already pass the angular cut.")

# ================================================================ K05: threshold scans
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
den = acc & has_bdpt
for thr, col in zip((0.05, 0.1, 0.2, 0.3, 0.5), C):
    eff_vs_pt(axs[0], kpt, np.abs(bdpt) < thr, den, K_EDGES, f"|Δpt/pt| < {thr}", col)
axs[0].set_title("momentum threshold scan, given dR < 0.02", fontsize=11)
den = acc & has_bdr & in_loose
for thr, col in zip((0.005, 0.01, 0.02, 0.05, 0.1), C):
    eff_vs_pt(axs[1], kpt, bdr < thr, den, K_EDGES, f"ΔR < {thr}", col)
axs[1].set_title("angular threshold scan, given |Δpt/pt| < 0.1", fontsize=11)
for ax in axs:
    style_pt_axis(ax, K_EDGES); ax.set_ylim(0.5, 1.02); ax.set_ylabel("efficiency"); ax.legend(fontsize=8.5, loc="lower right")
save(fig, "K05_kaon_threshold_scan.png",
     "Efficiency of alternative thresholds per gen pT bin. Left: momentum cut scanned for tracks already inside the angular cone. "
     "Right: angular cut scanned for momentum-compatible tracks. The producer uses 0.1 and 0.02.")

# ================================================================ K06/K07: 1D residuals per coarse pT bin
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
for (lo, hi), col in zip(COARSE, C):
    s = acc & has_bdpt & (kpt >= lo) & (kpt < hi)
    axs[0].hist(np.clip(bdpt[s], -0.5, 0.5), bins=np.linspace(-0.5, 0.5, 101), histtype="step", lw=1.6, color=col,
                density=True, label=f"{lo:g}–{hi:g} GeV  (n={s.sum()})")
    s = acc & has_bdr & in_loose & (kpt >= lo) & (kpt < hi)
    axs[1].hist(np.clip(bdr[s], 1e-4, 1), bins=np.geomspace(1e-4, 1, 80), histtype="step", lw=1.6, color=col, density=True,
                label=f"{lo:g}–{hi:g} GeV  (n={s.sum()})")
for y in (-PROD_DPT, PROD_DPT):
    axs[0].axvline(y, color=INK2, lw=1, ls="--")
axs[0].set_yscale("log"); axs[0].set_xlabel("(pt_reco − pt_gen)/pt_gen, best track in dR<0.02"); axs[0].set_ylabel("density")
axs[1].axvline(PROD_DR, color=INK2, lw=1, ls="--"); axs[1].set_xscale("log"); axs[1].set_yscale("log")
axs[1].set_xlabel("ΔR, best momentum-compatible track"); axs[1].set_ylabel("density")
for ax in axs:
    ax.legend(fontsize=8.5)
axs[0].set_title("Kaon momentum residual by gen pT", fontsize=11); axs[1].set_title("Kaon angular residual by gen pT", fontsize=11)
save(fig, "K06_kaon_residuals_1d_by_pt.png",
     "Normalised residual distributions per coarse gen-pT bin. Dashed: the producer's cuts. Tails outside the cut are the loss; "
     "compare their size between the softest and hardest bins.")

# ================================================================ K08: resolution summary
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
cx = centers(K_EDGES)
def q_in_bins(x, y, sel, q):
    return np.array([np.percentile(y[sel & (x >= lo) & (x < hi)], q) if np.sum(sel & (x >= lo) & (x < hi)) > 20 else np.nan
                     for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:])])
s = acc & has_bdpt
p16, p50, p84 = (q_in_bins(kpt, bdpt, s, q) for q in (16, 50, 84))
axs[0].errorbar(cx, 0.5 * (p84 - p16), xerr=xerr(K_EDGES), fmt="o", color=C[0], label="robust σ = (p84−p16)/2")
axs[0].errorbar(cx, np.abs(p50), xerr=xerr(K_EDGES), fmt="s", color=C[1], label="|median|")
axs[0].errorbar(cx, q_in_bins(kpt, np.abs(bdpt), s, 99), xerr=xerr(K_EDGES), fmt="^", color=C[2], label="p99 of |Δpt/pt|")
axs[0].axhline(PROD_DPT, color=INK2, ls="--", lw=1); axs[0].text(0.16, PROD_DPT * 1.1, "cut 0.1", color=INK2, fontsize=9)
axs[0].set_yscale("log"); axs[0].set_ylabel("|Δpt/pt|"); axs[0].set_title("Kaon momentum resolution (tracks in dR<0.02)", fontsize=11)
s = acc & has_bdr & in_loose
for q, col, mk in ((50, C[0], "o"), (84, C[1], "s"), (99, C[2], "^")):
    axs[1].errorbar(cx, q_in_bins(kpt, bdr, s, q), xerr=xerr(K_EDGES), fmt=mk, color=col, label=f"p{q}")
axs[1].axhline(PROD_DR, color=INK2, ls="--", lw=1); axs[1].text(0.16, PROD_DR * 1.1, "cut 0.02", color=INK2, fontsize=9)
axs[1].set_yscale("log"); axs[1].set_ylabel("ΔR"); axs[1].set_title("Kaon angular residual (pt-compatible tracks)", fontsize=11)
for ax in axs:
    style_pt_axis(ax, K_EDGES); ax.legend(fontsize=9, loc="upper right")
save(fig, "K08_kaon_resolution_summary.png",
     "Quantiles of the two residuals per gen pT bin against the producer's thresholds. A cut sits safely only when p99 is well below it.")

# ================================================================ K09: mis-association / ambiguity
fig, ax = plt.subplots(figsize=(8.5, 5))
den = acc & in_loose & has_trk
nearest_pass = near_dr_ok & near_dpt_ok
eff_vs_pt(ax, kpt, anypass & ~nearest_pass, den, K_EDGES, "a DIFFERENT track passes while the nearest fails (mis-association risk)", C[1])
eff_vs_pt(ax, kpt, K["n_pass"] > 1, den, K_EDGES, "more than one track passes both cuts (ambiguous)", C[3], "s")
eff_vs_pt(ax, kpt, ~anypass & in_loose, den, K_EDGES, "no track passes but one exists within 0.1 (resolution loss)", C[0], "^")
style_pt_axis(ax, K_EDGES); ax.set_yscale("log"); ax.set_ylim(1e-4, 1); ax.set_ylabel("fraction of gen kaons with a track within 0.1")
ax.legend(fontsize=8.5, loc="lower left"); ax.set_title("Resolution loss vs mis-association", fontsize=11)
save(fig, "K09_kaon_misassociation.png",
     "Separating the two failure mechanisms: a resolution loss (a track exists but fails the cuts) against mis-association "
     "(the nearest track fails but another one passes, or several pass).")

# ================================================================ MUONS
mpt = MU["gen_pt"]; macc = np.abs(MU["gen_eta"]) < 2.4
m_has = MU["muon_near_dr"] > SENT
m_loose = MU["muon_n_in_loose"] > 0
m_pass = MU["muon_any_pass"] == 1
m_bdpt = MU["muon_best_dpt_within_dr"]; m_has_bdpt = m_bdpt > SENT
m_bdr = MU["muon_best_dr_passing_dpt"]; m_has_bdr = m_bdr > SENT

fig, axs = plt.subplots(1, 2, figsize=(13, 5))
s = macc & m_has & (MU["muon_near_dr"] < PROD_DR)
h = axs[0].hist2d(mpt[s], MU["muon_near_dpt"][s], bins=[np.geomspace(1.5, 30, 50), np.linspace(-0.15, 0.15, 81)], cmap="Blues", norm=LogNorm())
quantile_curves(axs[0], mpt[s], MU["muon_near_dpt"][s], MU_EDGES)
for y in (-PROD_DPT, PROD_DPT):
    axs[0].axhline(y, color=C[1], lw=1.4)
style_pt_axis(axs[0], np.array([1.5, 30]), "gen muon $p_T$ [GeV]"); axs[0].set_ylabel("(pt_reco − pt_gen)/pt_gen, nearest Muon (dR<0.02)")
axs[0].set_title("Muon momentum residual vs gen pT", fontsize=11); axs[0].legend(fontsize=9, loc="lower right")
s = macc & m_has_bdr
h = axs[1].hist2d(mpt[s], np.clip(m_bdr[s], 1e-4, 2), bins=[np.geomspace(1.5, 30, 50), np.geomspace(1e-4, 2, 70)], cmap="Blues", norm=LogNorm())
axs[1].set_yscale("log"); quantile_curves(axs[1], mpt[s], np.clip(m_bdr[s], 1e-4, 2), MU_EDGES, qs=(50, 84, 97.5))
axs[1].axhline(PROD_DR, color=C[1], lw=1.4)
style_pt_axis(axs[1], np.array([1.5, 30]), "gen muon $p_T$ [GeV]"); axs[1].set_ylabel("ΔR, pt-compatible Muon"); axs[1].set_title("Muon angular residual vs gen pT", fontsize=11)
save(fig, "M01_muon_residuals_2d.png",
     "Same as K01/K02 for the two muon legs (μ+ and μ− combined), against the Muon table the J/ψ legs are built from.")

fig, ax = plt.subplots(figsize=(8.5, 5.2))
den = macc & m_loose
eff_vs_pt(ax, mpt, m_pass, den, MU_EDGES, "P(pass both | Muon within 0.1)", C[0])
eff_vs_pt(ax, mpt, m_has_bdr & (m_bdr < PROD_DR), den, MU_EDGES, "P(pt-compatible Muon has dR<0.02 | Muon)", C[1], "s")
eff_vs_pt(ax, mpt, m_has_bdpt & (np.abs(m_bdpt) < PROD_DPT), macc & m_has_bdpt, MU_EDGES, "P(|Δpt/pt|<0.1 | dR<0.02)", C[5], "D")
eff_vs_pt(ax, mpt, m_loose, macc, MU_EDGES, "P(a same-charge Muon within 0.1 | acceptance)", C[2], "^")
style_pt_axis(ax, MU_EDGES, "gen muon $p_T$ [GeV]"); ax.set_ylim(0, 1.05); ax.set_ylabel("efficiency"); ax.legend(fontsize=8.5, loc="lower right")
ax.set_title("Muon matching efficiency", fontsize=11)
save(fig, "M02_muon_eff.png", "Muon-leg matching efficiencies and their components, per gen muon pT.")

# ================================================================ CANDIDATE LEVEL
ckpt = CD["kaon_gen_pt"]; cacc = np.abs(CD["kaon_gen_eta"]) < 2.4
found = CD["found"] == 1
sig = found & (CD["prod_signal"] == 1)
gate = found & (CD["gate"] == 1)
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
eff_vs_pt(axs[0], ckpt, CD["have_tracks"] == 1, cacc, K_EDGES, "all three legs have a loose track", C[2], "^")
eff_vs_pt(axs[0], ckpt, found, cacc, K_EDGES, "a BuJpsiK candidate uses exactly those tracks", C[0])
eff_vs_pt(axs[0], ckpt, gate, cacc, K_EDGES, "… and passes the joint-CVH gate", C[3], "s")
axs[0].set_title("From gen decay to reconstructed truth candidate", fontsize=11); axs[0].set_ylabel("fraction of gen decays (|η_K|<2.4)")
eff_vs_pt(axs[1], ckpt, CD["prod_leg2"] == 1, found, K_EDGES, "producer matched the kaon leg", C[1], "s")
eff_vs_pt(axs[1], ckpt, CD["prod_nmatched"] == 3, found, K_EDGES, "producer matched all three legs", C[0])
eff_vs_pt(axs[1], ckpt, CD["prod_signal"] == 1, found, K_EDGES, "labelled exclusive signal (histmaker logic)", C[5], "D")
axs[1].set_title("Producer verdict on TRUE signal candidates", fontsize=11); axs[1].set_ylabel("fraction of truth candidates")
for ax in axs:
    style_pt_axis(ax, K_EDGES); ax.set_ylim(0, 1.05); ax.legend(fontsize=8.5, loc="lower right")
save(fig, "C01_candidate_level.png",
     "Left: how many gen decays end up as a reconstructed candidate built from the right three tracks. "
     "Right: of those true-signal candidates, the fraction the producer's gen matching labels as signal — the template's own efficiency vs kaon pT.")

# C02: the same vs RECO kaon pT (what the fit bins in)
fig, ax = plt.subplots(figsize=(8.5, 5.2))
rk = CD["reco_kpt"]
eff_vs_pt(ax, rk, CD["prod_signal"] == 1, found & gate & (rk > 0), K_EDGES, "labelled signal | truth candidate, gated", C[5], "D")
eff_vs_pt(ax, rk, CD["prod_leg2"] == 1, found & gate & (rk > 0), K_EDGES, "kaon leg matched | truth candidate, gated", C[1], "s")
style_pt_axis(ax, K_EDGES, "RECO kaon $p_T$ [GeV]"); ax.set_ylim(0.8, 1.02); ax.set_ylabel("fraction"); ax.legend(fontsize=9, loc="lower right")
ax.set_title("Template efficiency against the fit's own pT variable", fontsize=11)
save(fig, "C02_template_eff_vs_reco_pt.png",
     "The fraction of gated true-signal candidates that enter the signal template, binned in reconstructed kaon pT (the axis the fit uses).")

# C03: mass shapes, labelled-signal vs lost, per pT band, with ratio panels.
# The truth candidate's kaon track is the loose nearest one (dR<0.1, no pt cut), so a
# "lost" candidate is either the real kaon failing the producer's cuts or a wrong track
# picked up because the kaon was not reconstructed. Split on momentum compatibility.
plaus = (np.abs(CD["k_near_dpt"]) < 0.2) & (CD["k_near_dr"] < 0.05) & (CD["k_near_dr"] > SENT)
bands = [(0.15, 1.0), (1.0, 2.0), (2.0, 20.0)]
medges = np.linspace(5.0, 5.5, 51)
mc = 0.5 * (medges[1:] + medges[:-1])
fig, axs = plt.subplots(2, 3, figsize=(16, 7), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.05))
for j, (lo, hi) in enumerate(bands):
    s_all = found & gate & (ckpt >= lo) & (ckpt < hi)
    s_sig = s_all & (CD["prod_signal"] == 1)
    s_lp = s_all & (CD["prod_signal"] != 1) & plaus
    s_lw = s_all & (CD["prod_signal"] != 1) & ~plaus
    top, bot = axs[0, j], axs[1, j]
    hs, _ = np.histogram(CD["mass"][s_sig], medges); hp, _ = np.histogram(CD["mass"][s_lp], medges); hw, _ = np.histogram(CD["mass"][s_lw], medges)
    ns, npl, nw = hs.sum(), hp.sum(), hw.sum()
    top.errorbar(mc, hs / max(ns, 1), yerr=np.sqrt(hs) / max(ns, 1), fmt="o", ms=3.5, color=C[0], label=f"labelled signal (n={ns})")
    if npl > 0:
        top.errorbar(mc, hp / npl, yerr=np.sqrt(hp) / npl, fmt="s", ms=3.5, color=C[1], label=f"lost, kaon track plausible |Δpt/pt|<0.2 (n={npl})")
    if nw > 0:
        top.errorbar(mc, hw / nw, yerr=np.sqrt(hw) / nw, fmt="^", ms=3.5, color=C[2], label=f"lost, likely a different track (n={nw})")
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (hp / npl) / (hs / ns) if ns > 0 and npl > 0 else np.full_like(mc, np.nan)
        re = r * np.sqrt(np.where(hp > 0, 1 / hp, 0) + np.where(hs > 0, 1 / hs, 0))
    bot.errorbar(mc, r, yerr=re, fmt="o", ms=3, color=INK); bot.axhline(1, color=INK2, lw=0.8)
    bot.set_ylim(0, 4); bot.set_xlabel("m(μμK), joint CVH [GeV]"); bot.set_ylabel("plausible-lost / labelled")
    top.set_title(f"gen kaon pT {lo:g}–{hi:g} GeV", fontsize=11); top.set_yscale("log"); top.legend(fontsize=8)
axs[0, 0].set_ylabel("normalised")
save(fig, "C03_mass_labelled_vs_lost.png",
     "Gated joint-CVH mass of true-signal candidates: the producer-labelled template (blue), candidates the producer did not label whose kaon track is "
     "momentum-compatible with the gen kaon (orange: the real loss), and those whose kaon track is not (green: a different track picked up because the "
     "kaon was not reconstructed; correctly rejected). Ratio: orange over blue.")

# C04: the template as built vs full truth, inclusive and pT>1
fig, axs = plt.subplots(2, 2, figsize=(12, 7), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.05))
rows = []
for j, (lab, extra) in enumerate((("all kaon pT", np.ones_like(found, bool)), ("reco kaon pT > 1 GeV (fit selection)", rk > 1.0))):
    s_all = found & gate & extra
    s_fair = s_all & (plaus | (CD["prod_signal"] == 1))
    s_sig = s_all & (CD["prod_signal"] == 1)
    top, bot = axs[0, j], axs[1, j]
    ha, _ = np.histogram(CD["mass"][s_all], medges); hf, _ = np.histogram(CD["mass"][s_fair], medges); hs, _ = np.histogram(CD["mass"][s_sig], medges)
    top.errorbar(mc, ha, yerr=np.sqrt(ha), fmt="^", ms=3.5, color=C[2], label=f"all truth candidates, gated (n={ha.sum()})")
    top.errorbar(mc, hf, yerr=np.sqrt(hf), fmt="o", ms=3.5, color=C[0], label=f"… with a momentum-compatible kaon track (n={hf.sum()})")
    top.errorbar(mc, hs, yerr=np.sqrt(hs), fmt="s", ms=3.5, color=C[1], label=f"template as built (labelled) (n={hs.sum()})")
    with np.errstate(divide="ignore", invalid="ignore"):
        r = hs / hf; re = np.sqrt(r * (1 - r) / np.where(hf > 0, hf, 1))
        r2 = hs / ha; re2 = np.sqrt(r2 * (1 - r2) / np.where(ha > 0, ha, 1))
    bot.errorbar(mc, r2, yerr=re2, fmt="^", ms=3, color=C[2], label="template / all")
    bot.errorbar(mc, r, yerr=re, fmt="o", ms=3, color=INK, label="template / momentum-compatible")
    bot.axhline(1, color=INK2, lw=0.8); bot.set_ylim(0.6, 1.05); bot.legend(fontsize=8, loc="lower left")
    bot.set_xlabel("m(μμK), joint CVH [GeV]"); bot.set_ylabel("ratio")
    top.set_yscale("log"); top.set_title(lab, fontsize=11); top.legend(fontsize=8.5)
    for name, s_ in (("truth, all", s_all), ("truth, momentum-compatible kaon", s_fair), ("template", s_sig)):
        m = CD["mass"][s_]; q = np.percentile(m, [16, 50, 84]) if s_.sum() else [np.nan] * 3
        rows.append((lab, name, int(s_.sum()), q[1], 0.5 * (q[2] - q[0]), np.mean((m < 5.2) | (m > 5.36)) if s_.sum() else np.nan))
axs[0, 0].set_ylabel("candidates / 10 MeV")
save(fig, "C04_template_vs_truth.png",
     "The signal template as the histmaker builds it (orange) against all gated truth candidates (green) and against the fair reference: truth candidates "
     "whose kaon track is momentum-compatible with the gen kaon (blue). The black ratio is the template's genuine mass-dependent efficiency.")

# ================================================================ G: where the soft kaons go
jpsi_exists = (CD["mup_near_dr"] > SENT) & (CD["mup_near_dr"] < 0.1) & (CD["mum_near_dr"] > SENT) & (CD["mum_near_dr"] < 0.1)
k_track = (CD["k_near_dr"] > SENT) & (CD["k_near_dr"] < 0.1)
lab_sig = found & (CD["prod_signal"] == 1)
fig, (top, bot) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.3], hspace=0.05))
fine = np.geomspace(0.15, 20, 36)
tight = (CD["mup_tight"] == 1) & (CD["mum_tight"] == 1)
series = [("gen decays, |η_K|<2.4", cacc, C[5], "-"),
          ("… and both muons have a track (loose J/ψ exists)", cacc & jpsi_exists, C[0], "-"),
          ("… and both muons pass the tight proxy (global ∧ tracker ∧ nMatches>1) — the B+ channel's J/ψ", cacc & jpsi_exists & tight, INK2, "-"),
          ("… and the kaon has a track within ΔR<0.1 (entered a candidate)", cacc & jpsi_exists & tight & k_track, C[1], "-"),
          ("… and a BuJpsiK candidate uses the three tracks", cacc & jpsi_exists & tight & found, C[2], "-"),
          ("… and it passes the joint-CVH gate", cacc & jpsi_exists & tight & found & gate, C[3], "-"),
          ("… and the producer labels it signal (the template)", cacc & jpsi_exists & tight & lab_sig & gate, C[4], "-")]
ref, _ = np.histogram(ckpt[cacc & jpsi_exists & tight], fine)
for lab, m, col, ls in series:
    h, _ = np.histogram(ckpt[m], fine)
    top.step(fine[:-1], h, where="post", color=col, lw=1.6, label=lab)
    if m is not series[0][1]:
        with np.errstate(divide="ignore", invalid="ignore"):
            r = h / ref; re = np.sqrt(r * (1 - r) / np.where(ref > 0, ref, 1))
        bot.errorbar(np.sqrt(fine[:-1] * fine[1:]), r, yerr=re, fmt="o", ms=3, color=col)
top.set_yscale("log"); top.set_ylabel("decays / bin"); top.legend(fontsize=8.5, loc="lower right")
top.set_title("Gen kaon pT spectrum through the reconstruction chain", fontsize=11)
bot.set_ylabel("fraction of decays\nwith a TIGHT J/ψ"); bot.set_ylim(0, 1.05); bot.axvline(1.0, color=INK2, ls=":", lw=0.8)
bot.set_xscale("log"); bot.set_xlim(0.15, 20); bot.set_xlabel("gen kaon $p_T$ [GeV]")
save(fig, "G1_kaon_pt_cascade.png",
     "Where the soft kaons go. Each curve is the gen kaon pT spectrum after one more step of the chain; the ratio panel is relative to decays whose two muons "
     "pass the tight-selector proxy the B+ channel's J/ψ requires (so the muon loss is taken out of the kaon curves). The big step is 'the kaon has a track at all': the nano Track table holds only tracks that entered an "
     "AlCaReco candidate, so this is the AlCaReco B+ producer's bachelor selection (pT>0.1, |η|<2.5, m(μμK) in [5.0,5.5], DOCA(K, μ tracks) < 0.03 cm).")

# gate terms
fig, ax = plt.subplots(figsize=(9, 5.2))
okk = found & (CD["ok"] == 1) & (CD["mass"] > SENT)
with np.errstate(divide="ignore", invalid="ignore"):
    c2n = np.where(CD["ndof"] > 0, CD["chisq"] / np.where(CD["ndof"] > 0, CD["ndof"], 1), np.inf)
eff_vs_pt(ax, ckpt, found & (CD["ok"] != 1), found, K_EDGES, "jointCvhOk != 1 (refit did not converge)", C[0])
eff_vs_pt(ax, ckpt, okk & ~((CD["edm"] >= 0) & (CD["edm"] < 1e-2)), okk, K_EDGES, "edmRef ≥ 1e-2  | Ok", C[1], "s")
eff_vs_pt(ax, ckpt, okk & ~((CD["masserr"] > 0) & (CD["masserr"] < 1)), okk, K_EDGES, "massErr ∉ (0,1)  | Ok", C[2], "^")
eff_vs_pt(ax, ckpt, okk & ~((CD["chisq"] > 0) & (CD["chisq"] < 1e6) & (CD["ndof"] > 0) & (CD["ndof"] < 200)), okk, K_EDGES, "chisq<1e6 & ndof<200 (current gate)  | Ok", C[3], "v")
eff_vs_pt(ax, ckpt, okk & (c2n > 3), okk, K_EDGES, "chi2/ndof > 3  | Ok  (not a gate today)", C[5], "D")
eff_vs_pt(ax, ckpt, found & ~gate, found, K_EDGES, "fails the full gate", INK, "x")
style_pt_axis(ax, K_EDGES); ax.set_yscale("log"); ax.set_ylim(2e-4, 1); ax.set_ylabel("fraction of truth candidates lost")
ax.legend(fontsize=8.5, loc="upper right"); ax.set_title("Joint-CVH gate: which term removes true signal, per gen kaon pT", fontsize=11)
save(fig, "G2_gate_terms_vs_pt.png",
     "Loss of true-signal candidates to each term of the joint-CVH gate, per gen kaon pT. The chisq<1e6 and ndof<200 terms remove nothing; a chi2/ndof cut "
     "is shown for comparison although it is not applied today.")

# chi2/ndof distributions
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
for (lo, hi), col in zip(COARSE, C):
    s_ = okk & (ckpt >= lo) & (ckpt < hi)
    axs[0].hist(np.clip(c2n[s_], 0, 8), bins=np.linspace(0, 8, 81), histtype="step", lw=1.6, color=col, density=True, label=f"{lo:g}–{hi:g} GeV (n={s_.sum()})")
    axs[1].hist(np.clip(CD["ndof"][s_], 20, 90), bins=np.arange(20, 91, 1), histtype="step", lw=1.6, color=col, density=True, label=f"{lo:g}–{hi:g} GeV")
axs[0].set_yscale("log"); axs[0].set_xlabel("joint-CVH chi2 / ndof"); axs[0].set_ylabel("density"); axs[0].legend(fontsize=8.5)
axs[0].set_title("True-signal candidates: chi2/ndof by gen kaon pT", fontsize=11)
axs[1].set_xlabel("joint-CVH ndof"); axs[1].set_ylabel("density"); axs[1].legend(fontsize=8.5); axs[1].set_title("ndof by gen kaon pT (varies with hit count)", fontsize=11)
save(fig, "G3_chi2ndof_by_pt.png",
     "chi2/ndof and ndof of the joint CVH refit for true-signal candidates per gen kaon pT. ndof spans 30–70, so a cut on raw chisq is a cut on hit count in disguise; "
     "chi2/ndof is the quantity to gate on.")

# kaon track quality -- counts, not quantiles
kk = found & (CD["k_d0err"] > SENT)
ptb = np.geomspace(0.15, 20, 40)
fig, axs = plt.subplots(1, 3, figsize=(17.5, 5))
h = axs[0].hist2d(ckpt[kk], np.clip(CD["k_d0err"][kk], 1e-3, 0.3), bins=[ptb, np.geomspace(1e-3, 0.3, 60)], cmap="Blues", norm=LogNorm())
axs[0].set_yscale("log"); axs[0].axhline(0.03, color=C[1], lw=1.4); axs[0].text(2.2, 0.033, "AlCaReco DOCA(K, μ) cut = 0.03 cm", color=C[1], fontsize=9)
axs[0].set_ylabel("σ(d0) of the kaon track [cm]  (track-fit covariance)"); axs[0].set_title("Kaon impact-parameter uncertainty", fontsize=11)
fig.colorbar(h[3], ax=axs[0], label="candidates")
h = axs[1].hist2d(ckpt[kk], CD["k_hits"][kk], bins=[ptb, np.arange(2.5, 31.5, 1)], cmap="Blues", norm=LogNorm())
axs[1].set_ylabel("nValidHits of the kaon track"); axs[1].set_title("Kaon track length", fontsize=11)
fig.colorbar(h[3], ax=axs[1], label="candidates")
rel = np.where(CD["reco_kpt"] > 0, CD["k_pterr"] / np.where(CD["reco_kpt"] > 0, CD["reco_kpt"], 1), np.nan)
h = axs[2].hist2d(ckpt[kk], np.clip(rel[kk], 3e-3, 0.3), bins=[ptb, np.geomspace(3e-3, 0.3, 60)], cmap="Blues", norm=LogNorm())
axs[2].set_yscale("log"); axs[2].axhline(PROD_DPT, color=C[1], lw=1.4); axs[2].text(2.2, 0.11, "gen-match cut |Δpt/pt| < 0.1", color=C[1], fontsize=9)
axs[2].set_ylabel("σ(pt)/pt of the kaon track  (track-fit covariance)"); axs[2].set_title("Kaon momentum uncertainty", fontsize=11)
fig.colorbar(h[3], ax=axs[2], label="candidates")
for ax in axs:
    style_pt_axis(ax, K_EDGES)
save(fig, "G4_kaon_track_quality_vs_pt.png",
     "Track-fit quality of the kaon in true-signal candidates that survived to the nano, as counts. σ(d0) and σ(pt) are the track fit's own covariance "
     "(Track_d0Err, Track_ptErr), not residuals against truth. Left: the survivors' σ(d0) reaches the AlCaReco's 300 µm DOCA cut below ~0.7 GeV.")

# G6: the producer's DOCA cut, recomputed on the survivors
dk = found & (CD["doca_max"] > SENT)
fig, axs = plt.subplots(1, 2, figsize=(14, 5.2))
h = axs[0].hist2d(ckpt[dk], np.clip(CD["doca_max"][dk], 1e-4, 0.05), bins=[ptb, np.linspace(0, 0.035, 71)], cmap="Blues", norm=LogNorm())
axs[0].axhline(0.03, color=C[1], lw=1.4); axs[0].text(0.16, 0.0305, "cut: DOCA < 0.03 cm to each muon", color=C[1], fontsize=9)
style_pt_axis(axs[0], K_EDGES); axs[0].set_ylabel("max over the two muons of straight-line DOCA(K track, μ track) [cm]")
axs[0].set_title("Producer DOCA, recomputed for surviving true-signal candidates", fontsize=11); fig.colorbar(h[3], ax=axs[0], label="candidates")
for (lo, hi), col in zip(COARSE, C):
    s_ = dk & (ckpt >= lo) & (ckpt < hi)
    axs[1].hist(np.clip(CD["doca_max"][s_], 0, 0.035), bins=np.linspace(0, 0.035, 71), histtype="step", lw=1.6, color=col, density=True,
                label=f"{lo:g}–{hi:g} GeV (n={s_.sum()}, {100*np.mean(CD['doca_max'][s_] > 0.024):.0f}% in the top 20% below the cut)")
axs[1].axvline(0.03, color=INK2, ls="--", lw=1.2); axs[1].set_yscale("log"); axs[1].set_xlabel("max DOCA(K, μ) [cm]"); axs[1].set_ylabel("density")
axs[1].legend(fontsize=8.5); axs[1].set_title("DOCA shape per gen kaon pT: a distribution truncated at the cut is one the cut is eating", fontsize=10.5)
save(fig, "G6_producer_doca_survivors.png",
     "The AlCaReco producer's own selection variable (straight-line track-track DOCA, its formula rebuilt from Track_dxy/dz/phi/eta), for the true-signal "
     "candidates that survived it. Above ~2 GeV the distribution has died off long before 0.03 cm; below 1 GeV it is still flat at the cut, so the cut removes "
     "a large fraction there. The legend quotes the fraction of survivors in the last 20% of the allowed range.")

# G7: infer the DOCA cut's acceptance from the survivors' truncated half-Gaussian
def qb(y, sel, q):
    return np.array([np.percentile(y[sel & (ckpt >= lo) & (ckpt < hi)], q) if np.sum(sel & (ckpt >= lo) & (ckpt < hi)) > 20 else np.nan
                     for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:])])
cx = centers(K_EDGES)
from scipy.optimize import curve_fit
from scipy.special import erf
def max2_pdf(x, A, sig):
    """|d| = max of two independent half-Gaussians of width sig: pdf = 2 F f."""
    z = x / (sig * np.sqrt(2))
    return A * erf(z) * np.exp(-0.5 * (x / sig) ** 2)
def max2_acc(sig, cut=0.03):
    return erf(cut / (sig * np.sqrt(2))) ** 2
SIG_MAX = 0.5
db = np.linspace(0, 0.03, 31); dc = 0.5 * (db[1:] + db[:-1])
sig_fit, sig_err, acc_fit, acc_lo, acc_hi, at_bound = [], [], [], [], [], []
for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:]):
    s_ = dk & (ckpt >= lo) & (ckpt < hi)
    hh, _ = np.histogram(CD["doca_max"][s_], db)
    try:
        popt, pcov = curve_fit(max2_pdf, dc, hh, p0=[hh.max() if hh.max() > 0 else 1, 0.01], sigma=np.sqrt(np.maximum(hh, 1)),
                               bounds=([0, 1e-3], [np.inf, SIG_MAX]))
        sg, se = popt[1], np.sqrt(pcov[1, 1])
    except Exception:
        sg, se = np.nan, np.nan
    bound = np.isfinite(sg) and sg > 0.8 * SIG_MAX
    at_bound.append(bound)
    sig_fit.append(sg); sig_err.append(se if not bound else 0.0)
    acc_fit.append(max2_acc(sg) if np.isfinite(sg) else np.nan)
    acc_lo.append(max2_acc(sg + se) if np.isfinite(se) and not bound else 0.0)
    acc_hi.append(max2_acc(max(sg - se, 1e-3)) if np.isfinite(se) else np.nan)
sig_fit, acc_fit, acc_lo, acc_hi, at_bound = map(np.array, (sig_fit, acc_fit, acc_lo, acc_hi, at_bound))
fig, axs = plt.subplots(1, 2, figsize=(14, 5.2))
ok_ = ~at_bound
axs[0].errorbar(cx[ok_], sig_fit[ok_], yerr=np.array(sig_err)[ok_], xerr=xerr(K_EDGES)[:, ok_], fmt="o", color=C[0],
                label="σ_eff of the K–μ DOCA (fit: max of two half-Gaussians, d < 0.03)")
if at_bound.any():
    axs[0].errorbar(cx[at_bound], np.full(at_bound.sum(), 0.03), xerr=xerr(K_EDGES)[:, at_bound], yerr=[np.zeros(at_bound.sum()), np.full(at_bound.sum(), 0.03)],
                    fmt="o", color=C[0], mfc="white", lolims=True, label="distribution flat to the cut: σ_eff ≳ 0.03 (lower bound)")
axs[0].errorbar(cx, qb(CD["k_d0err"], kk, 50), xerr=xerr(K_EDGES), fmt="s", color=C[1], label="median σ(d0) of the surviving kaon track, for comparison")
axs[0].axhline(0.03, color=INK2, ls="--", lw=1); axs[0].text(2.2, 0.032, "cut 0.03 cm", color=INK2, fontsize=9)
axs[0].set_yscale("log"); axs[0].set_ylabel("[cm]"); axs[0].set_title("Width of the K–μ DOCA distribution vs gen kaon pT", fontsize=11)
axs[1].fill_between(cx[ok_], acc_lo[ok_], acc_hi[ok_], color=C[0], alpha=0.25)
axs[1].plot(cx[ok_], acc_fit[ok_], "o-", color=C[0], label="DOCA-cut acceptance implied by σ_eff: erf(0.03/σ√2)²")
if at_bound.any():
    axs[1].errorbar(cx[at_bound], np.full(at_bound.sum(), max2_acc(0.03)), xerr=xerr(K_EDGES)[:, at_bound], yerr=[np.full(at_bound.sum(), max2_acc(0.03)), np.zeros(at_bound.sum())],
                    fmt="o", color=C[0], mfc="white", uplims=True, label="σ_eff only bounded: acceptance ≲ this (upper limit)")
sd0 = qb(CD["k_d0err"], kk, 50)
axs[1].plot(cx, erf(0.03 / (sd0 * np.sqrt(2))), "D--", color=C[2], label="if the DOCA error ≈ surviving kaons' median σ(d0): erf(0.03/σ√2)  (1D, correlated limit)")
eff_vs_pt(axs[1], ckpt, k_track, cacc & jpsi_exists & tight, K_EDGES, "observed: kaon has a track, given a TIGHT J/ψ", C[1], "s")
eff_vs_pt(axs[1], ckpt, k_track, cacc & jpsi_exists, K_EDGES, "observed: kaon has a track, given a loose J/ψ", C[3], "^")
axs[1].set_ylim(0, 1.05); axs[1].set_ylabel("fraction"); axs[1].legend(fontsize=8, loc="center right", bbox_to_anchor=(1.0, 0.45))
axs[1].set_title("How much of the low-pT loss the DOCA cut explains", fontsize=11)
for ax in axs:
    style_pt_axis(ax, K_EDGES)
axs[0].legend(fontsize=8.5, loc="upper right")
save(fig, "G7_doca_acceptance_estimate.png",
     "Left: σ_eff from a fit to the survivors' max-DOCA below the cut with the correct model (maximum of two half-Gaussians, whose density rises linearly "
     "from zero). Where the distribution is still flat at the cut the fit only bounds σ from below (open markers). Right: the acceptance that σ_eff implies "
     "for the 0.03 cm cut, against the observed fraction of decays whose kaon entered any candidate. Where the blue sits above the orange, the gap is "
     "everything else — tracking inefficiency, kaon decay in flight, the m(μμK) window.")
rows_doca = list(zip(K_EDGES[:-1], K_EDGES[1:], sig_fit, acc_fit))

# reco kaon pT of the template
fig, (top, bot) = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.05))
h_all, _ = np.histogram(rk[found & gate], fine); h_lab, _ = np.histogram(rk[lab_sig & gate], fine)
top.step(fine[:-1], h_all, where="post", color=C[0], lw=1.6, label=f"true signal, gated (n={h_all.sum()})")
top.step(fine[:-1], h_lab, where="post", color=C[1], lw=1.6, label=f"labelled signal, gated = template (n={h_lab.sum()})")
top.axvline(1.0, color=INK2, ls="--", lw=1.2); top.text(1.05, h_all.max() * 0.6, "histmaker default --bachelorMinPt 1.0\n(since 2026-08-25)", fontsize=9, color=INK2)
below = (rk[lab_sig & gate] < 1.0).mean()
top.text(0.16, h_all.max() * 0.05, f"{100*below:.1f}% of the template sits below 1 GeV", fontsize=10, color=INK)
top.set_yscale("log"); top.set_ylabel("candidates / bin"); top.legend(fontsize=9, loc="lower right"); top.set_title("Reconstructed kaon pT of true-signal candidates in the nano", fontsize=11)
with np.errstate(divide="ignore", invalid="ignore"):
    r = h_lab / h_all; re = np.sqrt(r * (1 - r) / np.where(h_all > 0, h_all, 1))
bot.errorbar(np.sqrt(fine[:-1] * fine[1:]), r, yerr=re, fmt="o", ms=3, color=INK); bot.axhline(1, color=INK2, lw=0.8); bot.set_ylim(0.7, 1.02)
bot.set_xscale("log"); bot.set_xlim(0.15, 20); bot.set_xlabel("RECO kaon $p_T$ [GeV]"); bot.set_ylabel("template / truth")
save(fig, "G5_reco_kaon_pt_template.png",
     "The nano does contain true signal below 1 GeV of reco kaon pT; every histogram since the aem change has excluded it through the histmaker's "
     "--bachelorMinPt 1.0 default, not through the nano.")

# ================================================================ summary table
lines = [f"# Gen-matching study — summary ({nfiles} files, {nev} events, {nK} exclusive B+ -> J/psi(mumu) K decays)", "",
         "Producer thresholds: dR < 0.02 AND |pt_reco - pt_gen|/pt_gen < 0.1, same charge, status 1, smallest dR.", "",
         "## Kaon, per gen pT bin", "",
         "| pT [GeV] | N (|η|<2.4) | P(track in 0.1) | P(pass \\| track) | P(dR<0.02 \\| track) | P(dpt ok \\| dR ok) | robust σ(Δpt/pt) | p99 |Δpt/pt| | p99 ΔR (pt-ok, in 0.1) | misassoc. |",
         "|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:]):
    b = acc & (kpt >= lo) & (kpt < hi); n = b.sum()
    if n == 0:
        continue
    t = b & in_loose
    ptrk = t.sum() / n
    ppass = anypass[t].mean() if t.sum() else np.nan
    pdr = near_dr_ok[t].mean() if t.sum() else np.nan
    d2 = b & has_bdpt
    pdpt = (np.abs(bdpt[d2]) < PROD_DPT).mean() if d2.sum() else np.nan
    q = np.percentile(bdpt[d2], [16, 84]) if d2.sum() > 20 else [np.nan, np.nan]
    p99 = np.percentile(np.abs(bdpt[d2]), 99) if d2.sum() > 20 else np.nan
    d3 = b & has_bdr & in_loose
    p99r = np.percentile(bdr[d3], 99) if d3.sum() > 20 else np.nan
    mis = (anypass & ~(near_dr_ok & near_dpt_ok))[t].mean() if t.sum() else np.nan
    lines.append(f"| {lo:g}–{hi:g} | {n} | {ptrk:.3f} | {ppass:.3f} | {pdr:.3f} | {pdpt:.3f} | {0.5*(q[1]-q[0]):.4f} | {p99:.3f} | {p99r:.4f} | {mis:.4f} |")
lines += ["", "## Candidate level, per gen kaon pT bin", "",
          "| pT [GeV] | gen decays | truth cand. found | gated | producer: kaon leg matched | all 3 matched | labelled signal | lost & plausible kaon | lost & other track |",
          "|---|--:|--:|--:|--:|--:|--:|--:|--:|"]
for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:]):
    b = cacc & (ckpt >= lo) & (ckpt < hi); f = b & found
    if b.sum() == 0:
        continue
    lines.append(f"| {lo:g}–{hi:g} | {b.sum()} | {f.sum()} | {(f & gate).sum()} | "
                 f"{(CD['prod_leg2'][f]==1).mean() if f.sum() else np.nan:.4f} | {(CD['prod_nmatched'][f]==3).mean() if f.sum() else np.nan:.4f} | "
                 f"{(CD['prod_signal'][f]==1).mean() if f.sum() else np.nan:.4f} | "
                 f"{((CD['prod_signal'][f]!=1) & plaus[f]).mean() if f.sum() else np.nan:.4f} | {((CD['prod_signal'][f]!=1) & ~plaus[f]).mean() if f.sum() else np.nan:.4f} |")
lines += ["", "## Joint-CVH gate terms, loss of true-signal candidates per gen kaon pT bin", "",
          "| pT [GeV] | found | Ok!=1 | edmRef>=1e-2 (given Ok) | massErr out (given Ok) | chisq/ndof caps (given Ok) | chi2/ndof>3 (given Ok, not applied) | full gate |",
          "|---|--:|--:|--:|--:|--:|--:|--:|"]
for lo, hi in zip(K_EDGES[:-1], K_EDGES[1:]):
    b = cacc & (ckpt >= lo) & (ckpt < hi); f = b & found; o = f & (CD["ok"] == 1) & (CD["mass"] > SENT)
    if f.sum() == 0:
        continue
    def fr(num, den):
        return f"{num.sum()/den.sum():.4f}" if den.sum() else "nan"
    lines.append(f"| {lo:g}–{hi:g} | {f.sum()} | {fr(f & (CD['ok']!=1), f)} | {fr(o & ~((CD['edm']>=0)&(CD['edm']<1e-2)), o)} | "
                 f"{fr(o & ~((CD['masserr']>0)&(CD['masserr']<1)), o)} | {fr(o & ~((CD['chisq']>0)&(CD['chisq']<1e6)&(CD['ndof']>0)&(CD['ndof']<200)), o)} | "
                 f"{fr(o & (c2n>3), o)} | {fr(f & ~gate, f)} |")
lines += ["", "## DOCA cut: sigma_eff from the survivors and the implied acceptance", "", "| pT [GeV] | sigma_eff [cm] (max2 fit) | acc. erf(0.03/sigma sqrt2)^2 | acc. from median sigma(d0), 1D | observed P(K track \\| loose J/psi) | observed P(K track \\| tight J/psi) |", "|---|--:|--:|--:|--:|--:|"]
for (lo, hi, sg, a) in rows_doca:
    b = cacc & jpsi_exists & (ckpt >= lo) & (ckpt < hi)
    flag = " (bound)" if sg > 0.8 * SIG_MAX else ""
    bt = b & tight; i_ = int(np.searchsorted(K_EDGES, lo))
    lines.append(f"| {lo:g}–{hi:g} | {sg:.4f}{flag} | {a:.3f}{flag} | {erf(0.03/(sd0[i_]*np.sqrt(2))):.3f} | {k_track[b].mean() if b.sum() else np.nan:.3f} | {k_track[bt].mean() if bt.sum() else np.nan:.3f} |")
lines += ["", "## Template vs truth (gated joint-CVH mass)", "", "| selection | set | n | median [GeV] | (p84−p16)/2 [GeV] | fraction outside [5.20, 5.36] |", "|---|---|--:|--:|--:|--:|"]
for lab, name, n, med, w, tail in rows:
    lines.append(f"| {lab} | {name} | {n} | {med:.4f} | {w:.4f} | {tail:.4f} |")
summary = "\n".join(lines)
with open(os.path.join(OUT, "summary.md"), "w") as fh:
    fh.write(summary + "\n")
print(summary)

# ================================================================ index.html
html = ["<html><head><meta charset='utf-8'><title>gen-matching study</title>",
        "<style>body{font-family:sans-serif;max-width:1400px;margin:auto;color:#0b0b0b} img{max-width:100%;border:1px solid #ddd} "
        "figure{margin:24px 0} figcaption{color:#52514e;font-size:14px} pre{background:#f6f6f4;padding:12px;overflow-x:auto}</style></head><body>",
        f"<h1>Gen-matching study — B+ → J/ψ(μμ) K, {nfiles} MC files</h1>",
        "<p>Every exclusive decay found from the gen record; each leg confronted with reco objects with the producer's criteria "
        "(same charge, ΔR &lt; 0.02, |Δpt/pt| &lt; 0.1) left free. Unweighted counts.</p>"]
for name, cap in FIGS:
    html.append(f"<figure><img src='{name}'><figcaption><b>{name}</b> — {cap}</figcaption></figure>")
html.append("<h2>Summary table</h2><pre>" + summary.replace("<", "&lt;") + "</pre></body></html>")
with open(os.path.join(OUT, "index.html"), "w") as fh:
    fh.write("\n".join(html))
print("wrote index.html")
