"""Why the 2018 dimuon and muon cuts cannot flatten this combinatorial.

Each of those cuts was designed for a sample whose background was different. Here
the discriminating power is measured directly: signal against the no-b-ancestor
combinatorial, for every variable the 2018 selection cuts on. A cut can only help
if the two distributions differ, and mostly they do not -- the background contains
a REAL J/psi from a real b hadron, so every dimuon-level quantity looks like
signal by construction. What differs is the bachelor track, which is why the
useful handles are all on the kaon.

    python scripts/btojpsik/make_cvhnano_dimuon_figs.py <histmaker.hdf5>
"""

import argparse
import os

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt

from wremnants.postprocessing import rabbit_cvhnano_helpers as H
from wremnants.production.btojpsik_cvhnano_axes import GEN_NO_B_ANCESTOR, GEN_SIGNAL

MC = "InclusiveBToJpsiX_2016PostVFP"
SURFACE, INK_MUTED = "#fcfcfb", "#6b6b6b"
C_SIG, C_BKG = "#2a78d6", "#9c9c9c"

# (histogram, axis, label, 2018 cut value, cut direction)
PANELS = [
    ("nominal_cand_dimuonRawPt", "dimuon_raw_pt", r"$p_{T}(\mu\mu)$ [GeV]", 7.0, ">"),
    ("nominal_BuJpsiK_kvfCvhDimuonVtxProb", "kvfCvh_dimuon_vtxProb",
     "dimuon vertex probability", 0.1, ">"),
    ("nominal_BuJpsiK_kvfCvhDimuonAlphaBS", "kvfCvh_dimuon_alphaBS",
     r"dimuon $\alpha_{BS}$ [rad]", 0.4, "<"),
    ("nominal_cand_mu0Pt", "mu0_pt", r"$p_{T}(\mu_{0})$ [GeV]", 4.0, ">"),
]


def get(o, n):
    h = o[n]
    return h.get() if hasattr(h, "get") else h


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("infile")
    p.add_argument("-o", "--outdir",
                   default="/home/submit/pmlugato/public_html/slides/assets_cvhnano_truth")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    with H.open_results(args.infile) as r:
        o = r[MC]["output"]
        for hname, axname, label, cutval, direction in PANELS:
            if hname not in o:
                print(f"  absent, skipping: {hname}")
                continue
            h = get(o, hname)
            edges = h.axes[0].edges
            centers = h.axes[0].centers
            s = h.values()[:, GEN_SIGNAL].astype(float)
            b = h.values()[:, GEN_NO_B_ANCESTOR].astype(float)
            if s.sum() <= 0 or b.sum() <= 0:
                continue
            sn, bn = s / s.sum(), b / b.sum()
            keep = centers > cutval if direction == ">" else centers < cutval
            es, eb = sn[keep].sum(), bn[keep].sum()

            fig, ax = plt.subplots(figsize=(7.2, 4.2))
            fig.patch.set_facecolor(SURFACE)
            ax.set_facecolor(SURFACE)
            ax.stairs(bn, edges, color=C_BKG, lw=0, fill=True, alpha=0.85,
                      label="no b ancestor (combinatorial)")
            ax.stairs(sn, edges, color=C_SIG, lw=2.2, baseline=None,
                      label=r"signal B$^+\to$J/$\psi$K$^+$")
            ax.axvline(cutval, color="#c0392b", lw=1.8, ls="--")
            ax.annotate(f"2018 cut\n{direction} {cutval:g}", xy=(cutval, ax.get_ylim()[1]*0.86),
                        xytext=(6, 0), textcoords="offset points", fontsize=8.5,
                        color="#c0392b", va="top")
            ax.set_xlabel(label)
            ax.set_ylabel("fraction of candidates / bin")
            ax.set_xlim(edges[0], edges[-1])
            ax.set_ylim(0.0, None)
            ax.grid(alpha=0.22)
            ax.set_axisbelow(True)
            ax.legend(fontsize=8.5, frameon=False)
            ax.set_title(f"keeps {es:.1%} of signal, {eb:.1%} of combinatorial",
                         fontsize=9, color=INK_MUTED, loc="right")
            fig.tight_layout()
            out = f"fig_dimuon_{axname}.png"
            fig.savefig(f"{args.outdir}/{out}", bbox_inches="tight", dpi=140)
            plt.close(fig)
            print(f"wrote {out}  (signal {es:.3f}, combinatorial {eb:.3f})")


if __name__ == "__main__":
    main()
