"""Per-cell chi-square of a rabbit plot-ready file: data against the postfit
prediction, per (pT, eta, charge) cell over the mass bins, data variance only
(Neyman), so the sum over cells is comparable to the fit's saturated 2*deltaNLL.

rabbit's histogram plotter prints ONE chi-square per mapping (all cells at
once) on every per-cell figure; this gives each cell its own number and a
figure of chi2/ndf against eta for every (pT, charge) series.

    python scripts/btojpsik/cvhnano_cell_chi2.py <plotready.hdf5> -o <dir> [--fittype postfit|prefit]
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rabbit.io_tools as io

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
MUTED = "#6b6a63"
MAPPING = "Project btojpsik_cvhnano fit_bachelor_pt fit_bachelor_eta fit_bachelor_charge jointCvh_mass"
CHANNEL = "btojpsik_cvhnano"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("infile")
    p.add_argument("-o", "--outfolder", required=True)
    p.add_argument("--fittype", default="postfit", choices=["postfit", "prefit"])
    p.add_argument("--postfix", default="")
    p.add_argument("--noTitle", action="store_true")
    a = p.parse_args()
    os.makedirs(a.outfolder, exist_ok=True)

    fr = io.get_fitresult(a.infile)
    ch = fr["mappings"][MAPPING]["channels"][CHANNEL]
    keys = sorted(ch.keys())
    hd = ch["hist_data_obs"].get()
    hp = ch[f"hist_{a.fittype}_inclusive"].get()
    d, dv = np.asarray(hd.values(), float), np.asarray(hd.variances(), float)
    pr, pv = np.asarray(hp.values(), float), np.asarray(hp.variances(), float)
    if pv is None or not np.all(np.isfinite(pv)):
        pv = np.zeros_like(pr)
    axes = list(hp.axes.name)
    imass = axes.index("jointCvh_mass")
    # Data variance only (Neyman), so the sum over cells is comparable to the
    # fit's saturated 2*deltaNLL. Adding the postfit prediction band would
    # count the template statistics twice and hide the sparse cells.
    var = dv
    chi2 = np.where(var > 0, (d - pr) ** 2 / np.maximum(var, 1e-12), 0.0).sum(axis=imass)
    ndf = pr.shape[imass]
    pt_edges = hp.axes["fit_bachelor_pt"].edges
    eta_c = hp.axes["fit_bachelor_eta"].centers
    npt, neta, nq = chi2.shape
    lines = [f"# per-cell chi2/ndf, {a.fittype}, ndf = {ndf} mass bins per cell; keys {keys}",
             "| pT bin | eta bin | charge | chi2 | chi2/ndf |", "|---|---|---|---|---|"]
    for i in range(npt):
        for j in range(neta):
            for k in range(nq):
                lines.append(f"| {i} [{pt_edges[i]:.1f},{pt_edges[i+1]:.1f}] | {j} ({eta_c[j]:+.1f}) | {k} | {chi2[i,j,k]:.1f} | {chi2[i,j,k]/ndf:.2f} |")
    flat = chi2.ravel() / ndf
    lines.append(f"\nall cells: median chi2/ndf {np.median(flat):.2f}, mean {flat.mean():.2f}, max {flat.max():.2f}, "
                 f"cells above 2: {(flat > 2).sum()} of {flat.size}; total chi2 {chi2.sum():.1f} over {chi2.size * ndf} bins")
    txt = "\n".join(lines)
    with open(os.path.join(a.outfolder, f"cell_chi2_{a.fittype}{a.postfix}.md"), "w") as f:
        f.write(txt + "\n")
    print(txt.split("\n")[-1])

    fig, ax = plt.subplots(figsize=(12.0, 4.8))
    qlab = {0: "q = −", 1: "q = +"}
    for i in range(npt):
        for k in range(nq):
            col = PALETTE[(i * nq + k) % len(PALETTE)]
            ax.plot(eta_c + (k - 0.5) * 0.04, chi2[i, :, k] / ndf, "o-", ms=4, lw=1.4, color=col,
                    label=f"pT [{pt_edges[i]:.1f}, {pt_edges[i+1]:.1f}] GeV, {qlab.get(k, k)}")
    ax.axhline(1.0, color=MUTED, lw=0.8, ls="--")
    ax.set_xlabel("kaon eta (fit bin centre)")
    ax.set_ylabel(f"chi2 / {ndf} per cell")
    if not a.noTitle:
        ax.set_title(f"Per-cell chi2/ndf, {a.fittype}", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=8, ncol=3)
    ax.grid(True, color="#e6e5df", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    fn = os.path.join(a.outfolder, f"cell_chi2_{a.fittype}{a.postfix}.png")
    fig.savefig(fn, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(fn)
    return 0


if __name__ == "__main__":
    sys.exit(main())
