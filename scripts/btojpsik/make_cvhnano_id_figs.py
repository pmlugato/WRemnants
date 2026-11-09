"""Kaon-ID candidates: displacement and ionisation, in slices of bachelor pT.

The background this has to reject is now known: a real J/psi from a b hadron plus
a **soft prompt track from the primary vertex**, two thirds of them pions. Two
handles follow from that description, and both are studied against pT because
both depend on it.

**Displacement.** The signal kaon comes from the displaced B decay vertex; this
background comes from the PV. `Track_d0` referenced to the PV, and its
significance `d0/sigma(d0)`, measure exactly that difference. Neither existed
before this NanoAOD round (gaps G13 and G6).

**Ionisation.** dE/dx separates kaons from pions only below roughly 0.7 GeV/c,
which is where this background sits. It was dismissed earlier as "nearly
useless" on an inclusive distribution -- the wrong test, because inclusively the
sample is dominated by momenta where no separation is possible.

Signal is the simulated `BuJpsiK` category, whose bachelor is a true kaon by
construction. Background is `NoBAncestor`, which at low pT is 66.5% pions. So the
overlays are close to a kaon-against-pion comparison, and are labelled as what
they are rather than as truth-level particle identification.

    python scripts/btojpsik/make_cvhnano_id_figs.py <histmaker.hdf5>
"""

import argparse
import json
import os

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt

from wremnants.postprocessing import rabbit_cvhnano_helpers as H
from wremnants.production.btojpsik_cvhnano_axes import GEN_NO_B_ANCESTOR, GEN_SIGNAL

MC = "InclusiveBToJpsiX_2016PostVFP"
SURFACE, INK_MUTED = "#fcfcfb", "#6b6b6b"
# Four slices, chosen where the physics changes rather than evenly: 0.7 GeV is
# where ionisation stops separating K from pi, 1.0 is the 2018 production floor.
SLICES = [(0.0, 0.4), (0.4, 0.7), (0.7, 1.0), (1.5, 2.5)]
SLICE_COLOURS = ["#2a78d6", "#1baf7a", "#eb6834", "#4a3aa7"]

VARS = [
    ("quality_id_d0", "id_d0", r"$d_{0}(K)$ w.r.t. PV [cm]", "fig_id_d0.png", True),
    ("quality_id_d0Sig", "id_d0Sig", r"$d_{0}/\sigma_{d_{0}}$", "fig_id_d0sig.png", True),
    ("quality_id_dedx", "id_dedx", r"track d$E$/d$x$ (harmonic-2)", "fig_id_dedx.png", False),
    ("quality_id_dzPVSig", "id_dzPVSig", r"$d_{z}/\sigma_{d_{z}}$ w.r.t. PV",
     "fig_id_dzsig.png", True),
]


def get(o, n):
    h = o[n]
    return h.get() if hasattr(h, "get") else h


def slice_pt(h, axis, lo, hi, category):
    """Normalised distribution and its error, for one pT slice and category."""
    pt = h.axes["id_bachelor_pt"]
    m = (pt.centers >= lo) & (pt.centers < hi)
    iv = list(h.axes.name).index(axis)
    ipt = list(h.axes.name).index("id_bachelor_pt")
    ig = list(h.axes.name).index("genCategory")
    sl = [slice(None)] * 3
    sl[ipt] = m
    sl[ig] = category
    v = h.values()[tuple(sl)]
    var = h.variances()[tuple(sl)]
    ax_sum = tuple(i for i in range(v.ndim) if i != (0 if iv < ipt else v.ndim - 1))
    # After fancy-indexing the pT axis the variable axis keeps its relative order,
    # so sum everything that is not the variable axis rather than guessing.
    keep = 0 if iv < ipt else v.ndim - 1
    v = v.sum(axis=tuple(i for i in range(v.ndim) if i != keep))
    var = var.sum(axis=tuple(i for i in range(var.ndim) if i != keep))
    if v.sum() <= 0:
        return None, None, 0.0
    return v / v.sum(), np.sqrt(var) / v.sum(), float(v.sum())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("infile")
    p.add_argument("-o", "--outdir",
                   default="/home/submit/pmlugato/public_html/slides/assets_cvhnano_truth")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    summary = {}

    with H.open_results(args.infile) as r:
        o = r[MC]["output"]
        for hname, axname, label, out, signed in VARS:
            if hname not in o:
                print(f"  absent, skipping: {hname}")
                continue
            h = get(o, hname)
            edges = h.axes[axname].edges
            fig, ax = plt.subplots(figsize=(7.6, 4.5))
            fig.patch.set_facecolor(SURFACE)
            ax.set_facecolor(SURFACE)
            # Colour carries the pT slice and line style carries the category, so
            # the legend is 4 + 2 entries instead of 8. Eight entries on a
            # distribution plot is not a legend, it is a wall.
            handles = []
            for (lo, hi), col in zip(SLICES, SLICE_COLOURS):
                sv, se, sn = slice_pt(h, axname, lo, hi, GEN_SIGNAL)
                bv, be, bn = slice_pt(h, axname, lo, hi, GEN_NO_B_ANCESTOR)
                if sv is None or bv is None:
                    continue
                ax.stairs(sv, edges, color=col, lw=2.0, baseline=None)
                ax.stairs(bv, edges, color=col, lw=1.4, ls="--", baseline=None)
                handles.append(mpl.lines.Line2D([], [], color=col, lw=2.0,
                                                label=f"{lo:g}–{hi:g} GeV"))
                summary.setdefault(axname, {})[f"{lo}-{hi}"] = dict(
                    signal_n=sn, background_n=bn)
            handles += [
                mpl.lines.Line2D([], [], color=INK_MUTED, lw=2.0, label="signal"),
                mpl.lines.Line2D([], [], color=INK_MUTED, lw=1.4, ls="--",
                                 label="combinatorial"),
            ]
            ax.set_xlabel(label)
            ax.set_ylabel("fraction of candidates / bin")
            ax.set_xlim(edges[0], edges[-1])
            ax.set_ylim(0.0, None)
            ax.grid(alpha=0.22)
            ax.set_axisbelow(True)
            ax.legend(handles=handles, fontsize=8.5, frameon=False, ncol=2,
                      title="bachelor $p_{T}$", title_fontsize=8.5)
            fig.tight_layout()
            fig.savefig(f"{args.outdir}/{out}", bbox_inches="tight", dpi=140)
            plt.close(fig)
            print(f"wrote {out}")

        # --- what a cut on each variable could buy, per pT slice --------------
        #
        # Two cut families, because the two handles are different shapes. A
        # displacement cut keeps |x| ABOVE a threshold: the signal is displaced
        # and the background is not. An ionisation cut keeps x inside a BAND
        # around the kaon expectation. Both are scanned and the better reported.
        def roc(centers, sv, bv, kind):
            if kind == "absmin":
                ts = np.unique(np.abs(centers))
                keeps = [np.abs(centers) >= t for t in ts]
            else:
                mode = centers[int(np.argmax(sv))]
                ts = np.linspace(0.0, float(np.abs(centers - mode).max()), 300)
                keeps = [np.abs(centers - mode) <= t for t in ts]
            eff = np.array([sv[k].sum() for k in keeps])
            rej = np.array([1.0 - bv[k].sum() for k in keeps])
            o = np.argsort(eff)
            return eff[o], rej[o]

        def rej_at(eff, rej, target):
            if eff.max() < target:
                return float("nan")
            return float(np.interp(target, eff, rej))

        print(f"\n{'variable':12s}{'pT slice':>11s}{'family':>9s}"
              f"{'rej @ 90% eff':>15s}{'rej @ 80% eff':>15s}{'signal n':>10s}")
        print("-" * 72)
        rocs = {}
        for hname, axname, label, out, signed in VARS:
            if hname not in o:
                continue
            h = get(o, hname)
            centers = h.axes[axname].centers
            for lo, hi in SLICES:
                sv, _, sn = slice_pt(h, axname, lo, hi, GEN_SIGNAL)
                bv, _, bn = slice_pt(h, axname, lo, hi, GEN_NO_B_ANCESTOR)
                if sv is None or bv is None or sn < 20:
                    continue
                rows = {}
                for kind in ("absmin", "band"):
                    e, rj = roc(centers, sv, bv, kind)
                    rows[kind] = (rej_at(e, rj, 0.90), rej_at(e, rj, 0.80))
                pick = max(rows, key=lambda k: (rows[k][0] if rows[k][0] == rows[k][0] else -1))
                r90, r80 = rows[pick]
                print(f"{axname:12s}{f'{lo:g}-{hi:g}':>11s}{pick:>9s}"
                      f"{r90:15.3f}{r80:15.3f}{sn:10.0f}")
                rocs.setdefault(axname, {})[f"{lo}-{hi}"] = dict(
                    rej_at_90=r90, rej_at_80=r80, family=pick, signal_n=sn)
        summary["rejection"] = rocs

    json.dump(summary, open(f"{args.outdir}/id_summary.json", "w"), indent=1)
    print(f"\nwrote {args.outdir}/id_summary.json")


if __name__ == "__main__":
    main()
