"""Template variations of the B+ tensor, drawn as var/nominal against mass.

Reads the tensor through rabbit's FitDebugData (the same object
rabbit_debug_inputdata.py uses), so what is plotted is exactly what the fit
sees: the log-normal up and down templates of each systematic on the signal
process, in the chosen (pT, eta, charge) cells.

    python scripts/btojpsik/cvhnano_variation_plots.py <tensor.hdf5> -o <dir> \\
        --systs scaleA12 scalee12 scaleM12 sigRes_eta12 \\
        --cells 0,12,0 2,12,0 0,0,0 2,23,1

One figure per systematic family: a row of panels, one per cell, each panel
showing up/nominal and down/nominal with the nominal template's statistical
band. Ratios rather than absolute templates, because the variations are
percent-level shifts on a peak that spans two orders of magnitude.
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rabbit import debugdata, inputdata

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK = "#1a1a19"
MUTED = "#6b6a63"
CH = "btojpsik_cvhnano"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tensor")
    p.add_argument("-o", "--outfolder", required=True)
    p.add_argument("--process", default="BuJpsiK")
    p.add_argument("--systs", nargs="+", required=True, help="systematic names")
    p.add_argument("--cells", nargs="+", required=True, help="pt,eta,charge bin indices")
    p.add_argument("--postfix", default="")
    p.add_argument("--noTitle", action="store_true")
    p.add_argument("--ylim", type=float, default=0.3, help="symmetric ratio range about 1")
    p.add_argument("--minEvents", type=float, default=3.0,
                   help="bins whose nominal holds fewer simulated events than this (before the "
                   "global scale) are hatched and their ratios not drawn: a ratio of a single "
                   "scaled event is noise, not a variation")
    p.add_argument("--theta", type=float, default=1.0,
                   help="draw the variation at +-theta instead of +-1, using rabbit's own "
                   "interpolation (var/nominal)**theta; small theta shows what the fit sees "
                   "near its best value")
    p.add_argument("--signalScale", type=float, default=18.2433,
                   help="global template scale, to convert bin content back to simulated events")
    a = p.parse_args()
    os.makedirs(a.outfolder, exist_ok=True)

    ind = inputdata.FitInputData(a.tensor)
    dd = debugdata.FitDebugData(ind)
    nom = dd.nominal_hists[CH][{"processes": a.process}]
    sy = dd.syst_hists[CH]
    axes = list(nom.axes.name)
    mass = nom.axes["jointCvh_mass"]
    x = mass.centers
    edges = mass.edges
    pt_edges = nom.axes["fit_bachelor_pt"].edges
    eta_edges = nom.axes["fit_bachelor_eta"].edges
    cells = [tuple(int(v) for v in c.split(",")) for c in a.cells]

    def cell_label(c):
        i, j, k = c
        return (f"pT [{pt_edges[i]:.1f}, {pt_edges[i+1]:.1f}], "
                f"η [{eta_edges[j]:+.1f}, {eta_edges[j+1]:+.1f}], q {'−' if k == 0 else '+'}")

    def sel(h, c):
        i, j, k = c
        return h[{"fit_bachelor_pt": i, "fit_bachelor_eta": j, "fit_bachelor_charge": k}]

    out = []
    for s in a.systs:
        fig, axs = plt.subplots(1, len(cells), figsize=(3.3 * len(cells), 3.6), sharey=True)
        axs = np.atleast_1d(axs)
        for ax, c in zip(axs, cells):
            n = sel(nom, c)
            nv = np.asarray(n.values(), float)
            nvar = np.asarray(n.variances(), float)
            up = np.asarray(sel(sy[{"processes": a.process, "systs": s, "DownUp": "Up"}], c).values(), float)
            dn = np.asarray(sel(sy[{"processes": a.process, "systs": s, "DownUp": "Down"}], c).values(), float)
            ok = nv > a.minEvents * a.signalScale
            if a.theta != 1.0:  # log-normal interpolation, as in rabbit
                pos = (nv > 0) & (up > 0) & (dn > 0)
                up = np.where(pos, nv * (np.where(pos, up, 1) / np.where(pos, nv, 1)) ** a.theta, up)
                dn = np.where(pos, nv * (np.where(pos, dn, 1) / np.where(pos, nv, 1)) ** a.theta, dn)
            ru = np.where(ok, up / np.maximum(nv, 1e-12), np.nan)
            rd = np.where(ok, dn / np.maximum(nv, 1e-12), np.nan)
            rel = np.where(ok, np.sqrt(nvar) / np.maximum(nv, 1e-12), np.nan)
            ax.stairs(np.where(ok, 1 + rel, 1), edges, baseline=np.where(ok, 1 - rel, 1), fill=True,
                      color="#e6e5df", label="template stat.")
            sparse = ~ok
            if sparse.any():
                ax.stairs(np.where(sparse, 1 + a.ylim, 1), edges, baseline=np.where(sparse, 1 - a.ylim, 1),
                          fill=True, color="none", hatch="////", edgecolor="#c8c7bf", lw=0,
                          label=f"< {a.minEvents:g} MC events")
            # mean mass shift of the template under the +1 sigma variation, in MeV
            if nv.sum() > 0 and up.sum() > 0:
                dm = 1e3 * ((up * x).sum() / up.sum() - (nv * x).sum() / nv.sum())
                ax.text(0.98, 0.04, f"Δ⟨m⟩(+{a.theta:g}σ) = {dm:+.2f} MeV", transform=ax.transAxes,
                        ha="right", va="bottom", fontsize=8, color=MUTED)
            ax.stairs(ru, edges, color=PALETTE[0], lw=1.6, label=f"+{a.theta:g}σ / nominal")
            ax.stairs(rd, edges, color=PALETTE[1], lw=1.6, label=f"−{a.theta:g}σ / nominal")
            ax.axhline(1.0, color=MUTED, lw=0.8, ls="--")
            ax.set_title(cell_label(c), fontsize=9, color=INK, loc="left")
            ax.set_xlabel("m(μμK) [GeV]")
            ax.grid(True, color="#e6e5df", lw=0.6)
            ax.spines[["top", "right"]].set_visible(False)
            if a.ylim:
                ax.set_ylim(1 - a.ylim, 1 + a.ylim)
        axs[0].set_ylabel(f"{s}: variation / nominal")
        axs[0].legend(frameon=False, fontsize=8, loc="upper left")
        if not a.noTitle:
            fig.suptitle(f"{a.process}, systematic {s}", x=0.01, ha="left", fontsize=11, color=INK)
        fn = os.path.join(a.outfolder, f"variation_{s}{a.postfix}.png")
        fig.savefig(fn, dpi=150, bbox_inches="tight")
        plt.close(fig)
        out.append(fn)
        print(fn)
    return 0


if __name__ == "__main__":
    sys.exit(main())
