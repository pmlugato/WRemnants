"""Every variable the CVH NanoAOD B+ histmaker cuts on, with the cut drawn on it.

Reads a histmaker HDF5 produced with `--allaxes` and, for each selection in
`scripts/histmakers/btojpsik_cvhnano.py`, plots the distribution of the variable
that selection tests, data against stacked simulation with a ratio panel, and
marks the threshold.

The input should be the *loosened* run -- gates wide open, pT floors at zero --
because a variable plotted after its own cut is truncated exactly where the
interesting part is. The arm-ok flag and the sentinel rejection stay on even
there: a sentinel (-99, not NaN) passes a one-sided bound, so admitting them
would corrupt every distribution on the page rather than only its own.

Cuts that are declared but OFF at the nominal working point are plotted too, and
labelled as such -- they are the ones worth arguing about.

Example:

    python scripts/btojpsik/cvhnano_cut_variable_plots.py \
        --infile .../btojpsik_cvhnano_cutvars_loose.hdf5 \
        -o ~/public_html/mz/cutvars/
"""

import argparse
import os
import sys

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plotting"))

from btojpsik_cvhnano_datamc import (  # noqa: E402
    COLOR_OTHER_B,
    COLOR_SIGNAL,
    COLOR_UNMATCHED,
    DATA_SAMPLE,
    GEN_DATA,
    GEN_NO_B_ANCESTOR,
    GEN_REAL_B_BACKGROUND,
    GEN_SIGNAL,
    INK,
    INK_MUTED,
    MC_SAMPLE,
    SURFACE,
    load,
    project,
    style,
)

from wums import logging  # noqa: E402

logger = logging.child_logger(__name__)

CAND = "BuJpsiK"

# The colour reserved for the cut itself. Deliberately not one of the three
# categorical series colours -- the line is annotation, not data.
COLOR_CUT = "#b4232a"


class Cut:
    """One selection: which axis it tests, where it sits, whether it is on.

    `keep` is "above" or "below" for a one-sided bound, or None when the
    variable is cut on in a way a threshold cannot express (a flag, a pairing).
    """

    def __init__(
        self,
        var,
        label,
        keep=None,
        value=None,
        active=False,
        flag=None,
        logy=False,
        logx=False,
        note="",
    ):
        self.var = var
        self.label = label
        self.keep = keep
        self.value = value
        self.active = active
        self.flag = flag
        self.logy = logy
        self.logx = logx
        self.note = note


# The selection list, in the order `_selection_list` applies it. Values are the
# defaults recorded in the meta_info of the eta24joint production, which is the
# histogram every A/e/M result so far was built from.
CUTS = [
    # -- convergence gate ---------------------------------------------------
    Cut(f"{CAND}_jointCvhEdmRef", "jointCvh edmRef", "below", 1e-2, True,
        "--gateEdmRef", logy=True, logx=True,
        note="The only gate that bites: 15.7% of flagged-good candidates in data."),
    Cut(f"{CAND}_jointCvhChisq", "jointCvh chisq", "below", 1e6, True,
        "--gateChisq", logy=True,
        note="Threshold is far off the axis; this guards the pathological tail only."),
    Cut(f"{CAND}_jointCvhNdof", "jointCvh ndof", "below", 200.0, True,
        "--gateNdof", logy=True,
        note="Threshold is off the axis; an upper bound on an essentially fixed quantity."),
    Cut(f"{CAND}_jointCvhMassErr", "jointCvh massErr [GeV]", "below", 1.0, True,
        "--gateMassErr", logy=True,
        note="Threshold is off the axis; guards the pathological tail."),
    # -- muon legs ----------------------------------------------------------
    Cut("cand_mu0Pt", "leading muon pT [GeV]", "above", 3.0, True, "--muonMinPt"),
    Cut("cand_mu1Pt", "subleading muon pT [GeV]", "above", 3.0, True, "--muonMinPt"),
    Cut("cand_mu0Eta", "leading muon eta", None, 2.4, True, "--muonMaxEta",
        note="Two-sided |eta| < 2.4; the AlCaReco is already inside it."),
    Cut("cand_mu1Eta", "subleading muon eta", None, 2.4, True, "--muonMaxEta",
        note="Two-sided |eta| < 2.4; the AlCaReco is already inside it."),
    # -- bachelor -----------------------------------------------------------
    Cut(f"{CAND}_kaonPt", "bachelor pT [GeV]", "above", 1.0, True, "--bachelorMinPt",
        logy=True,
        note="The load-bearing cut: purity 0.214 -> 0.579 for 96.8% of the signal kept."),
    Cut(f"{CAND}_kaonEta", "bachelor eta", None, 2.4, True, "--bachelorMaxEta",
        note="Two-sided |eta| < 2.4."),
    # -- candidate ----------------------------------------------------------
    Cut(f"{CAND}_jointCvhMass", "jointCvh B mass [GeV]", None, None, True,
        "--massLow / --massHigh",
        note="Window [5.0, 5.5], set at AlCaReco production; the fit uses a subset."),
    Cut(f"{CAND}_nLegsRefit", "nLegsRefit", "above", 3, True, "--minLegsRefit",
        note="All three legs must carry a CVH refit."),
    # -- declared but OFF ---------------------------------------------------
    Cut(f"{CAND}_kaonPt", "bachelor pT [GeV], upper bound", "below", 8.0, False,
        "--bachelorMaxPt", logy=True,
        note="OFF. The 2018 analysis capped the kaon at 8 GeV."),
    Cut("cand_bachelorD0", "bachelor d0 [cm]", None, None, False,
        "--minBachelorAbsD0", logy=True,
        note="OFF. Two-sided |d0|; referenced to the origin, not the PV (gap G13)."),
    Cut("cand_bachelorD0Sig", "bachelor d0 / sigma", None, None, False,
        "--minBachelorAbsD0Sig", logy=True, note="OFF. Two-sided."),
    Cut("cand_bachelorDedxHarmonic2", "bachelor dE/dx (harmonic-2)", "above", None,
        False, "--minBachelorDedx",
        note="OFF. Separates K from pi only below ~0.7 GeV/c."),
    Cut("cand_bachelorNValidHits", "bachelor nValidHits", "above", None, False,
        "--bachelorMinValidHits"),
    Cut("cand_bachelorNormChi2", "bachelor track chi2/ndof", "below", None, False,
        "--maxChisqNdof", note="OFF. Track-level, distinct from the joint-fit chisq."),
    Cut(f"{CAND}_jointCvhNiter", "jointCvh nIter", "below", None, False,
        "--maxJointCvhNiter", logy=True),
    Cut(f"{CAND}_jointCvhSl3d", "jointCvh sl3d", "above", None, False, "--minSl3d",
        logy=True, note="OFF. Displacement significance, 3D."),
    Cut(f"{CAND}_jointCvhAlphaBS", "jointCvh alphaBS [rad]", "below", None, False,
        "--maxAlphaBS", logy=True,
        note="OFF. Transverse pointing angle; the 3D one is the right variable (gap G8)."),
    Cut(f"{CAND}_kvfCvhDimuonVtxProb", "dimuon vertex probability", "above", None,
        False, "--minDimuonVtxProb", logy=True,
        note="OFF. Borrowed from the kvfCvh arm: the joint arm has no dimuon vertex."),
    Cut(f"{CAND}_kvfCvhDimuonSl3d", "dimuon sl3d", "above", None, False,
        "--minDimuonSl3d", logy=True, note="OFF. Borrowed from the kvfCvh arm."),
    Cut(f"{CAND}_kvfCvhDimuonAlphaBS", "dimuon alphaBS [rad]", "below", None, False,
        "--maxDimuonAlphaBS", logy=True, note="OFF. Borrowed from the kvfCvh arm."),
    Cut(f"{CAND}_jpsiPt", "dimuon pT [GeV]", "above", None, False, "--minDimuonPt"),
    Cut(f"{CAND}_vertexChi2", "candidate vertex chi2", "below", None, False,
        "--minVtxProb", logy=True,
        note="OFF, and never applied in any run: the 2018 selection's tightest cut "
             "was the equivalent probability at 0.3."),
]


def cut_plot(data_hist, mc_hist, cut, outpath):
    """Data vs stacked simulation for one cut variable, with the cut marked."""
    edges, d_val, d_var = project(data_hist, [GEN_DATA])
    _, sig_val, sig_var = project(mc_hist, [GEN_SIGNAL])
    _, oth_val, oth_var = project(mc_hist, GEN_REAL_B_BACKGROUND)
    _, unm_val, unm_var = project(mc_hist, [GEN_NO_B_ANCESTOR])

    mc_val = sig_val + oth_val + unm_val
    mc_var = sig_var + oth_var + unm_var
    if mc_val.sum() <= 0 or d_val.sum() <= 0:
        logger.warning("Empty histogram for %s, skipping", cut.label)
        return None

    # Shape comparison only: the inclusive sample's cross section is a
    # placeholder, so the simulation is normalised to the data yield.
    scale = d_val.sum() / mc_val.sum()
    sig_val, oth_val, unm_val = sig_val * scale, oth_val * scale, unm_val * scale
    mc_val, mc_var = mc_val * scale, mc_var * scale**2

    centres = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    fig, (ax, rax) = plt.subplots(
        2, 1, figsize=(7.2, 6.0), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.07},
    )

    stack, labels, colors = [], [], []
    for vals, name, color in (
        (sig_val, r"MC: $\mathrm{B}^{+}\to\mathrm{J}/\psi\,\mathrm{K}^{+}$", COLOR_SIGNAL),
        (oth_val, "MC: other real b decay", COLOR_OTHER_B),
        (unm_val, "MC: no b ancestor", COLOR_UNMATCHED),
    ):
        if vals.sum() > 0:
            stack.append(vals)
            labels.append(name)
            colors.append(color)

    ax.hist([centres] * len(stack), bins=edges, weights=stack, stacked=True,
            color=colors, label=labels, histtype="stepfilled",
            edgecolor=SURFACE, linewidth=1.0)
    ax.errorbar(centres, d_val, yerr=np.sqrt(d_var), xerr=widths / 2, fmt="o",
                markersize=3.5, color=INK, ecolor=INK, elinewidth=1.0,
                capsize=0, linestyle="none", label="Data")

    # The cut itself. Drawn on both panels so the ratio can be read on the
    # kept side alone, and the rejected side shaded rather than only lined --
    # a bare line does not say which way the cut goes.
    drawn = False
    if cut.value is not None and cut.keep in ("above", "below"):
        lo, hi = ax.get_xlim()
        if lo < cut.value < hi:
            for a in (ax, rax):
                a.axvline(cut.value, color=COLOR_CUT, linewidth=1.4,
                          linestyle="--", zorder=5)
                a.axvspan(lo if cut.keep == "above" else cut.value,
                          cut.value if cut.keep == "above" else hi,
                          color=COLOR_CUT, alpha=0.07, zorder=0)
            ax.set_xlim(lo, hi)
            drawn = True
    elif cut.value is not None and cut.keep is None:
        # Two-sided |x| < value.
        lo, hi = ax.get_xlim()
        for a in (ax, rax):
            for v in (-cut.value, cut.value):
                if lo < v < hi:
                    a.axvline(v, color=COLOR_CUT, linewidth=1.4,
                              linestyle="--", zorder=5)
        ax.set_xlim(lo, hi)
        drawn = lo < cut.value < hi or lo < -cut.value < hi

    if cut.logy:
        ax.set_yscale("log")
        pos = np.concatenate([d_val[d_val > 0], mc_val[mc_val > 0]])
        if pos.size:
            ax.set_ylim(max(0.3, pos.min() * 0.5), pos.max() * 4)
    else:
        ax.set_ylim(bottom=0)
    if cut.logx:
        ax.set_xscale("log")
        rax.set_xscale("log")

    ax.set_ylabel("Candidates / bin")
    ax.grid(axis="y", alpha=0.7)
    ax.set_axisbelow(True)
    ax.legend(loc="best", fontsize=9)

    state = "ACTIVE" if cut.active else "declared, OFF"
    thr = ""
    if cut.value is not None:
        sign = {"above": ">", "below": "<"}.get(cut.keep, "|x| <")
        thr = f"  ·  {sign} {cut.value:g}"
        if not drawn:
            thr += "  (off axis)"
    ax.set_title(
        f"{cut.flag}  —  {state}{thr}", fontsize=10, color=INK, loc="left"
    )
    ax.text(0.0, 1.11,
            "B$^{+}\\to$J/$\\psi$K$^{+}$, CVH NanoAOD, 2016postVFP  —  loosened "
            "selection, MC normalised to data",
            transform=ax.transAxes, fontsize=8.5, color=INK_MUTED)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(mc_val > 0, d_val / mc_val, np.nan)
        rerr = np.where(mc_val > 0,
                        np.sqrt(d_var / mc_val**2 + (d_val**2) * mc_var / mc_val**4),
                        np.nan)
    rax.axhline(1.0, color=INK_MUTED, linewidth=0.9)
    rax.errorbar(centres, ratio, yerr=rerr, xerr=widths / 2, fmt="o",
                 markersize=3.5, color=INK, ecolor=INK, elinewidth=1.0,
                 capsize=0, linestyle="none")
    finite = np.isfinite(ratio)
    if finite.any():
        span = np.nanpercentile(np.abs(ratio[finite] - 1.0), 95)
        rax.set_ylim(1 - max(0.15, min(1.5, span * 1.6)),
                     1 + max(0.15, min(1.5, span * 1.6)))
    rax.set_ylabel("Data / MC", fontsize=9)
    rax.set_xlabel(cut.label)
    rax.grid(axis="y", alpha=0.7)
    rax.set_axisbelow(True)
    if cut.note:
        rax.text(0.0, -0.42, cut.note, transform=rax.transAxes,
                 fontsize=8.5, color=INK_MUTED, va="top", wrap=True)

    fig.savefig(outpath, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    return outpath


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--infile", required=True, help="Loosened histmaker HDF5")
    p.add_argument("-o", "--outfolder", required=True)
    args = p.parse_args()

    style()
    os.makedirs(args.outfolder, exist_ok=True)
    samples = load(args.infile)
    d, m = samples[DATA_SAMPLE], samples[MC_SAMPLE]

    written, missing = [], []
    for i, cut in enumerate(CUTS):
        key = f"nominal_{cut.var}"
        if key not in d or key not in m:
            missing.append((cut.flag, cut.var))
            continue
        tag = cut.flag.strip("-").replace(" / ", "_").replace("-", "_")
        name = f"cut_{i:02d}_{'on' if cut.active else 'off'}_{tag}.png"
        out = cut_plot(d[key], m[key], cut, os.path.join(args.outfolder, name))
        if out:
            written.append((out, cut))

    print(f"\nWrote {len(written)} plots to {args.outfolder}")
    for out, cut in written:
        print(f"  {'ON ' if cut.active else 'off'}  {os.path.basename(out)}")
    if missing:
        print(f"\n{len(missing)} cut variables had no histogram in the input:")
        for flag, var in missing:
            print(f"  {flag:28s} {var}")
    return written, missing


if __name__ == "__main__":
    main()
