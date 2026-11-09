#!/usr/bin/env python
"""Summary figures for the joint mZ + J/psi fit.

The central rabbit scripts cover prefit/postfit distributions
(`rabbit_plot_hists.py`), pulls and grouped impacts
(`rabbit_plot_pulls_and_impacts.py`), and the parameter covariance block
(`rabbit_plot_cov.py`). This adds the three views specific to the joint fit that
none of them produce:

1. how much each shared scale parameter is constrained, against eta
2. how strongly mZ correlates with each of them, against eta
3. what the grouped impacts on mZ are, as a chart rather than a table

Written for slides: ~2.5:1 aspect, no in-figure titles (the slide title carries
the wording), and nothing below readable size once scaled to a 11.5 in box.
"""

import argparse
import os

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt

from rabbit import io_tools

PARAMS = ["A", "e", "M"]
LABEL = {"A": "A  (b-field)", "e": "e  (material)", "M": "M  (alignment)"}
INK = "#1F1F1F"
GREY = "#595959"
# One hue per family, distinguishable in grey-scale by marker as well.
COL = {"A": "#3B6FB6", "e": "#C1553B", "M": "#4C8C5A"}
MARK = {"A": "o", "e": "s", "M": "^"}


def style(ax):
    ax.tick_params(labelsize=11, colors=INK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GREY)
    ax.grid(axis="y", color="#DDDDDD", lw=0.6)
    ax.set_axisbelow(True)


def load(path):
    fr = io_tools.get_fitresult(path)
    p, p0 = fr["parms"].get(), fr["parms_prefit"].get()
    names = [n.decode() if isinstance(n, bytes) else str(n) for n in p.axes["parms"]]
    cov = fr["cov"].get()
    cnames = [n.decode() if isinstance(n, bytes) else str(n) for n in cov.axes[0]]
    return (
        fr,
        {n: i for i, n in enumerate(names)},
        p.values(),
        np.sqrt(p.variances()),
        np.sqrt(p0.variances()),
        cov.values(),
        {n: i for i, n in enumerate(cnames)},
    )


def n_eta(idx):
    return sum(1 for k in idx if k.startswith("A") and k[1:].isdigit())


def fig_constraints(idx, err, err0, neta, out):
    fig, ax = plt.subplots(figsize=(11.0, 4.2))
    x = np.arange(neta)
    for fam in PARAMS:
        r = np.array([err[idx[f"{fam}{k}"]] / err0[idx[f"{fam}{k}"]] for k in range(neta)])
        ax.plot(x, r, MARK[fam] + "-", color=COL[fam], ms=5, lw=1.4, label=LABEL[fam])
    ax.axhline(1.0, color=GREY, ls="--", lw=1.2)
    ax.text(
        neta - 0.5, 1.005, "prefit width — data said nothing",
        ha="right", va="bottom", fontsize=10, color=GREY,
    )
    ax.set_xlabel("nuisance eta bin", fontsize=12, color=INK)
    ax.set_ylabel("postfit / prefit width", fontsize=12, color=INK)
    ax.set_ylim(0, 1.12)
    ax.set_xlim(-0.5, neta - 0.5)
    ax.legend(fontsize=11, frameon=False, ncol=3, loc="lower center")
    style(ax)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def fig_mz_corr(idx, cidx, cov, neta, out):
    fig, ax = plt.subplots(figsize=(11.0, 4.2))
    m = cidx["massShiftZ100MeV"]
    x = np.arange(neta)
    for fam in PARAMS:
        c = []
        for k in range(neta):
            j = cidx[f"{fam}{k}"]
            d = np.sqrt(cov[m, m] * cov[j, j])
            c.append(cov[m, j] / d if d > 0 else 0.0)
        ax.plot(x, c, MARK[fam] + "-", color=COL[fam], ms=5, lw=1.4, label=LABEL[fam])
    ax.axhline(0.0, color=GREY, ls="--", lw=1.2)
    ax.set_xlabel("nuisance eta bin", fontsize=12, color=INK)
    ax.set_ylabel(r"correlation with $m_Z$", fontsize=12, color=INK)
    ax.set_xlim(-0.5, neta - 0.5)
    ax.legend(fontsize=11, frameon=False, ncol=3, loc="upper center")
    style(ax)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def fig_impacts(fr, out, top=12):
    gi = fr["impacts_grouped"].get()
    gnames = [n.decode() if isinstance(n, bytes) else str(n) for n in gi.axes[-1]]
    pn = [n.decode() if isinstance(n, bytes) else str(n) for n in gi.axes[0]]
    row = pn.index("massShiftZ100MeV")
    vals = np.abs(gi.values()[row]) * 100.0  # NOI units -> MeV
    order = np.argsort(vals)[::-1][:top][::-1]

    fig, ax = plt.subplots(figsize=(10.0, 4.6))
    y = np.arange(len(order))
    cols = ["#3B6FB6" if gnames[o] == "scale_AeM" else "#9AA4AE" for o in order]
    ax.barh(y, vals[order], color=cols, height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels([gnames[o] for o in order], fontsize=11, color=INK)
    for yi, o in zip(y, order):
        ax.text(vals[o], yi, f"  {vals[o]:.2f}", va="center", fontsize=10, color=INK)
    ax.set_xlabel(r"impact on $m_Z$  [MeV]", fontsize=12, color=INK)
    ax.set_xlim(0, vals[order].max() * 1.18)
    style(ax)
    ax.grid(axis="x", color="#DDDDDD", lw=0.6)
    ax.grid(axis="y", visible=False)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def fig_pair_corr(cidx, cov, neta, out):
    """Adjacent-pair correlation: the signature of a too-fine nuisance basis."""
    fig, ax = plt.subplots(figsize=(11.0, 4.0))
    for fam in PARAMS:
        ks, cs = [], []
        for k in range(0, neta - 1, 2):
            a, b = cidx[f"{fam}{k}"], cidx[f"{fam}{k+1}"]
            d = np.sqrt(cov[a, a] * cov[b, b])
            if d > 0:
                ks.append(k)
                cs.append(cov[a, b] / d)
        ax.plot(ks, cs, MARK[fam] + "-", color=COL[fam], ms=5, lw=1.4, label=LABEL[fam])
    ax.axhline(0.0, color=GREY, ls="--", lw=1.2)
    ax.set_xlabel("first bin of the adjacent pair", fontsize=12, color=INK)
    ax.set_ylabel(r"$\rho$(bin $k$, bin $k{+}1$)", fontsize=12, color=INK)
    ax.set_ylim(-1.05, 1.05)
    ax.legend(fontsize=11, frameon=False, ncol=3, loc="upper center")
    style(ax)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fitresult")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    sfx = f"_{args.tag}" if args.tag else ""

    fr, idx, val, err, err0, cov, cidx = load(args.fitresult)
    neta = n_eta(idx)
    print(f"nuisance eta bins: {neta}  ({3*neta} shared scale parameters)")

    i = idx["massShiftZ100MeV"]
    print(f"massShiftZ100MeV = {val[i]:+.5f} +/- {err[i]:.5f}  -> sigma(mZ) = {err[i]*100:.2f} MeV")

    o = args.outdir
    fig_constraints(idx, err, err0, neta, f"{o}/joint_constraints{sfx}.png")
    fig_mz_corr(idx, cidx, cov, neta, f"{o}/joint_mz_correlation{sfx}.png")
    fig_impacts(fr, f"{o}/joint_impacts{sfx}.png")
    fig_pair_corr(cidx, cov, neta, f"{o}/joint_pair_correlation{sfx}.png")


if __name__ == "__main__":
    main()
