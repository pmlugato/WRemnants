#!/usr/bin/env python3
"""Plots for the AlCaReco -> CVH refit -> NanoAOD stream.

Reads the NanoAOD written by
`Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py ... nanoOut=<path>` and makes the
figures used in the status deck.

Palette: MIT-deck orange #E87722 + blue #5A7BB8. That pair is validated
(chroma floor, CVD separation dE 21.9 protan, normal-vision dE 27.8); the
orange sits below 3:1 against white, so every bar carries a direct value label
rather than relying on the fill alone.

Usage:
    plot_alcareco_nano.py <nano.root> <outdir>
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot

ORANGE = "#E87722"
BLUE = "#5A7BB8"
INK = "#111111"
MUTED = "#6B6B6B"
RULE = "#E5E5E5"

# PDG masses [GeV]
M_JPSI = 3.0969
M_BPLUS = 5.27934


def _style(ax):
    """Recessive axes: no top/right spine, hairline grid behind the marks."""
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(RULE)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(True, color=RULE, linewidth=0.8, alpha=0.9)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def _save(fig, outdir, name):
    path = os.path.join(outdir, name)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", path)


def mass_hist(vals, pdg, pdg_label, xlabel, title, outdir, name, xlim=None,
              show_offset=True):
    """Single-series distribution: one hue, no legend -- the title names it."""
    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    _style(ax)
    ax.hist(vals, bins=60, range=xlim, color=ORANGE, edgecolor="white", linewidth=0.4)
    ax.axvline(pdg, color=INK, linestyle="--", linewidth=1.4)
    ax.annotate(f"PDG {pdg_label}\n{pdg:.4f}", xy=(pdg, ax.get_ylim()[1] * 0.99),
                xytext=(5, 0), textcoords="offset points",
                color=INK, fontsize=9, ha="left", va="top")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("candidates")
    ax.set_title(title, color=INK, fontsize=11, loc="left", pad=8)
    med = np.median(vals)
    # one row per candidate: a J/psi shared by N bachelors enters N times
    stat = f"n = {len(vals)} cand.\nmedian = {med:.4f}"
    if show_offset:  # only meaningful where the distribution actually peaks
        stat += f"\noffset = {(med-pdg)*1000:+.1f} MeV"
    ax.annotate(stat,
                xy=(0.02, 0.95), xycoords="axes fraction",
                color=MUTED, fontsize=9, va="top")
    _save(fig, outdir, name)


def bars(labels, values, colors, ylabel, title, outdir, name, fmt="{:.1f}%"):
    """Categorical bars with direct value labels (the contrast relief)."""
    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    _style(ax)
    x = np.arange(len(labels))
    # 2px-equivalent gap between adjacent fills
    ax.bar(x, values, width=0.62, color=colors, edgecolor="white", linewidth=2)
    for xi, v in zip(x, values):
        ax.annotate(fmt.format(v), xy=(xi, v), xytext=(0, 4),
                    textcoords="offset points", ha="center",
                    color=INK, fontsize=10, fontweight="600")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, color=INK, fontsize=9)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, max(values) * 1.22)
    ax.set_title(title, color=INK, fontsize=11, loc="left", pad=8)
    _save(fig, outdir, name)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    infile, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    t = uproot.open(infile)["Events"]
    a = t.arrays(["BuJpsiK_mass", "BuJpsiK_corMass", "BuJpsiK_corMassErr",
                  "BuJpsiK_kaonPt", "nBuJpsiK", "nTrack", "nMuon"], library="np")
    flat = lambda k: np.concatenate([x for x in a[k] if len(x)])

    raw = flat("BuJpsiK_mass")
    cor = flat("BuJpsiK_corMass")
    kpt = flat("BuJpsiK_kaonPt")
    ok = cor > 0

    print(f"events={t.num_entries} candidates={len(raw)} corrected={ok.sum()}")
    print(f"raw B+  median={np.median(raw):.4f}")
    print(f"corMass median={np.median(cor[ok]):.5f}  (PDG J/psi {M_JPSI})")

    # 1. the refit dimuon mass -- the headline closure test
    mass_hist(cor[ok], M_JPSI, "J/$\\psi$",
              r"CVH-refit dimuon mass  $m(\mu\mu)$  [GeV]",
              "CVH-refit dimuon mass vs PDG J/$\\psi$",
              outdir, "cvh_dimuon_mass.png", xlim=(2.9, 3.3))

    # 2. raw B+ mass straight off the persisted candidate
    mass_hist(raw, M_BPLUS, "B$^+$", r"raw candidate mass  $m(\mu\mu K)$  [GeV]",
              "Raw B$^+$ mass from the AlCaReco candidate",
              outdir, "raw_bplus_mass.png", xlim=(5.0, 5.6), show_offset=False)

    # 3. bachelor kaon pT -- why plimit matters
    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    _style(ax)
    ax.hist(kpt, bins=50, range=(0, 4), color=ORANGE, edgecolor="white", linewidth=0.4)
    ax.axvline(1.0, color=INK, linestyle="--", linewidth=1.4)
    ax.annotate("old plimit\n1.0 GeV", xy=(1.0, ax.get_ylim()[1] * 0.92),
                xytext=(6, 0), textcoords="offset points",
                color=INK, fontsize=9, va="top")
    ax.set_xlabel(r"bachelor kaon $p_T$ [GeV]")
    ax.set_ylabel("candidates")
    ax.set_title("Bachelor kaons are soft -- most sat below the old floor",
                 color=INK, fontsize=11, loc="left", pad=8)
    frac = (kpt < 1.0).mean() * 100
    ax.annotate(f"{frac:.0f}% below 1.0 GeV", xy=(0.40, 0.72), xycoords="axes fraction",
                color=MUTED, fontsize=9)
    _save(fig, outdir, "kaon_pt.png")

    # 4. plimit effect on the single-track maker. Matched comparison: the SAME
    #    600 events / 1003 attempted bachelor tracks, only plimit changed.
    bars(["plimit = 1.0 GeV\n(old default)", "plimit = 0.05 GeV\n(new default)"],
         [428 / 1003 * 100, 936 / 1003 * 100], [BLUE, ORANGE],
         "single-track fits succeeded [%]",
         "Propagation floor: 1003 bachelor tracks, same events",
         outdir, "plimit_effect.png")

    # 5. convergence audit on the REFERENCE-BLOCK criterion (measured elsewhere)
    bars(["dimuon\nat cap", "single-track\nat cap", "single-track\nnever converge"],
         [1.3, 15.8, 8.2], [ORANGE, BLUE, BLUE],
         "fraction of successful fits [%]",
         "Reference-block EDM audit (cap = 10 iterations)",
         outdir, "convergence_audit.png")


if __name__ == "__main__":
    main()
