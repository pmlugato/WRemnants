"""Data/MC comparison and fit-arm overlay plots for the CVH NanoAOD histmaker.

Reads the HDF5 written by ``scripts/histmakers/btojpsik_cvhnano.py``. Every
comparison carries error bars and a ratio panel; a bare overlay of two
histograms hides exactly the disagreements these plots exist to find.

The simulation is normalised to the data yield. It has to be: the inclusive
b -> J/psi X sample carries a placeholder cross section, so only shapes are
being compared here, and the plots say so.

Example:

    python scripts/plotting/btojpsik_cvhnano_datamc.py \
        --infile btojpsik_cvhnano_v1.hdf5 -o plots/
"""

import argparse
import os

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from wums import ioutils, logging

mpl.use("Agg")

logger = logging.child_logger(__name__)

DATA_SAMPLE = "Charmonium_2016PostVFP_cvhNano"
MC_SAMPLE = "InclusiveBToJpsiX_2016PostVFP"

# genCategory encoding. Imported, never re-declared: this file used to mirror
# the integers by hand, and the encoding has since grown from four bins to
# eight, which a mirrored copy would have got silently wrong.
from wremnants.production.btojpsik_cvhnano_axes import (  # noqa: E402
    GEN_DATA,
    GEN_MC_CATEGORIES,
    GEN_NO_B_ANCESTOR,
    GEN_SIGNAL,
    GEN_TEMPLATE_CATEGORIES,
)

# Real b-hadron decays other than the signal. Small -- 0.46% of the gated peak
# window against 12.38% signal -- and each peaks 100-200 MeV below the B mass.
GEN_REAL_B_BACKGROUND = [c for c in GEN_TEMPLATE_CATEGORIES if c != GEN_SIGNAL]

# Categorical slots 1, 2, 3 of the reference palette, validated for CVD
# separation and normal-vision separation on a light surface.
COLOR_SIGNAL = "#2a78d6"
COLOR_OTHER_B = "#eb6834"
COLOR_UNMATCHED = "#1baf7a"
# Slots 1, 2, 7 for the three FIT arms: passes all-pairs CVD and
# normal-vision separation with every slot clear of 3:1 contrast.
# The raw pre-fit mass is not a fit arm, so it is drawn as a muted dashed
# reference rather than competing for a categorical hue.
ARM_COLORS = ["#2a78d6", "#eb6834", "#4a3aa7"]

# Ink, never a series colour.
INK = "#1c1c1c"
INK_MUTED = "#6b6b6b"
SURFACE = "#fcfcfb"

# Variables worth a data/MC comparison, and the axis label to print.
DEFAULT_VARIABLES = {
    "BuJpsiK_jointCvhMass": r"$m(\mu\mu K)$, joint CVH [GeV]",
    "BuJpsiK_kvfCvhMass": r"$m(\mu\mu K)$, KVF on CVH tracks [GeV]",
    "BuJpsiK_jointCvhPt": r"$p_{T}(B)$, joint CVH [GeV]",
    "BuJpsiK_kaonPt": r"$p_{T}(K)$ [GeV]",
    "BuJpsiK_kaonEta": r"$\eta(K)$",
    "cand_mu0Pt": r"$p_{T}(\mu_{0})$ [GeV]",
    "cand_mu0Eta": r"$\eta(\mu_{0})$",
    "cand_mu1Pt": r"$p_{T}(\mu_{1})$ [GeV]",
    # The joint arm's J/psi mass is a constrained parameter, identical for
    # every candidate, so it is deliberately not plotted. This is the raw
    # dimuon rebuilt from the muon tracks -- the only unconstrained one.
    "cand_dimuonRawMass": r"$m(\mu\mu)$, raw tracks [GeV]",
    "cand_dimuonRawPt": r"$p_{T}(\mu\mu)$, raw tracks [GeV]",
    "BuJpsiK_jointCvhJpsiPt": r"$p_{T}(J/\psi)$, joint CVH [GeV]",
    "BuJpsiK_kvfCvhVtxProb": r"$B$ vertex probability, KVF CVH",
    "BuJpsiK_kvfCvhDimuonVtxProb": r"dimuon vertex probability, KVF CVH",
    "BuJpsiK_kvfCvhSl3d": r"$L_{3D}/\sigma$ ($B$), KVF CVH",
    "BuJpsiK_kvfCvhAlphaBS": r"$\alpha_{BS}$ ($B$), KVF CVH",
    "BuJpsiK_jointCvhEdmRef": r"joint CVH edmRef",
    "nPV": r"$N_{PV}$",
    "nBuJpsiK": r"$N$ candidates per event",
    "cand_bachelorDedxHarmonic2": r"$K$ track d$E$/d$x$ (harmonic-2)",
}

# The three fit arms, for the overlay plot. The raw mass is handled separately.
RAW_MASS_VAR = ("BuJpsiK_mass", "raw (pre-fit)")
ARM_MASS_VARS = [
    ("BuJpsiK_kvfRawMass", "KVF, raw tracks"),
    ("BuJpsiK_kvfCvhMass", "KVF, CVH tracks"),
    ("BuJpsiK_jointCvhMass", "joint N-body CVH"),
]


def style():
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "axes.edgecolor": INK_MUTED,
            "axes.labelcolor": INK,
            "axes.linewidth": 0.8,
            "text.color": INK,
            "xtick.color": INK_MUTED,
            "ytick.color": INK_MUTED,
            "xtick.labelcolor": INK,
            "ytick.labelcolor": INK,
            "grid.color": "#e3e3e1",
            "grid.linewidth": 0.6,
            "font.size": 11,
            "legend.frameon": False,
        }
    )


def load(infile):
    """Read every histogram eagerly.

    The stored objects are lazy proxies onto the open HDF5 file, so they have
    to be resolved before the file is closed rather than carried out of scope.
    """
    out = {}
    with h5py.File(infile, "r") as handle:
        for sample in (DATA_SAMPLE, MC_SAMPLE):
            if sample not in handle:
                raise KeyError(
                    f"{sample} not in {infile}. Present: {list(handle.keys())}"
                )
            stored = ioutils.pickle_load_h5py(handle[sample])["output"]
            out[sample] = {name: stored[name].get() for name in stored.keys()}
    return out


def project(hist_obj, gen_bins):
    """Sum a histogram over the requested genCategory bins.

    Returns (edges, values, variances). Values are summed over the second axis,
    which is always genCategory in this histmaker's output.
    """
    values = hist_obj.values()
    variances = hist_obj.variances()
    if variances is None:
        variances = values
    sel = list(gen_bins)
    return (
        hist_obj.axes[0].edges,
        values[:, sel].sum(axis=1),
        variances[:, sel].sum(axis=1),
    )


def data_mc_plot(data_hist, mc_hist, label, outpath):
    edges, d_val, d_var = project(data_hist, [GEN_DATA])
    _, sig_val, sig_var = project(mc_hist, [GEN_SIGNAL])
    _, oth_val, oth_var = project(mc_hist, GEN_REAL_B_BACKGROUND)
    _, unm_val, unm_var = project(mc_hist, [GEN_NO_B_ANCESTOR])

    mc_val = sig_val + oth_val + unm_val
    mc_var = sig_var + oth_var + unm_var
    if mc_val.sum() <= 0 or d_val.sum() <= 0:
        logger.warning("Empty histogram for %s, skipping", label)
        return None

    # Shape comparison only: the simulation's cross section is a placeholder.
    scale = d_val.sum() / mc_val.sum()
    sig_val, oth_val, unm_val = sig_val * scale, oth_val * scale, unm_val * scale
    mc_val, mc_var = mc_val * scale, mc_var * scale**2

    centres = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.0, 6.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.07},
    )

    stack, labels, colors = [], [], []
    for vals, name, color in (
        (
            sig_val,
            r"MC: $\mathrm{B}^{+}\to\mathrm{J}/\psi\,\mathrm{K}^{+}$",
            COLOR_SIGNAL,
        ),
        (oth_val, "MC: other real b decay", COLOR_OTHER_B),
        (unm_val, "MC: no b ancestor", COLOR_UNMATCHED),
    ):
        if vals.sum() > 0:
            stack.append(vals)
            labels.append(name)
            colors.append(color)

    # A 2px surface-coloured edge separates adjacent stacked fills.
    ax.hist(
        [centres] * len(stack),
        bins=edges,
        weights=stack,
        stacked=True,
        color=colors,
        label=labels,
        histtype="stepfilled",
        edgecolor=SURFACE,
        linewidth=1.0,
    )
    ax.errorbar(
        centres,
        d_val,
        yerr=np.sqrt(d_var),
        xerr=widths / 2,
        fmt="o",
        markersize=3.5,
        color=INK,
        ecolor=INK,
        elinewidth=1.0,
        capsize=0,
        linestyle="none",
        label="Data",
    )
    ax.set_ylabel("Candidates / bin")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.7)
    ax.set_axisbelow(True)
    ax.legend(loc="best", fontsize=9)
    ax.set_title(
        "B$^{+}\\to$J/$\\psi$K$^{+}$, CVH NanoAOD, 2016postVFP  "
        "— MC normalised to data",
        fontsize=10,
        color=INK,
        loc="left",
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(mc_val > 0, d_val / mc_val, np.nan)
        # Data and simulation statistical errors, propagated onto the ratio.
        rerr = np.where(
            mc_val > 0,
            np.sqrt(d_var / mc_val**2 + (d_val**2) * mc_var / mc_val**4),
            np.nan,
        )
    rax.axhline(1.0, color=INK_MUTED, linewidth=0.9)
    rax.errorbar(
        centres,
        ratio,
        yerr=rerr,
        xerr=widths / 2,
        fmt="o",
        markersize=3.5,
        color=INK,
        ecolor=INK,
        elinewidth=1.0,
        capsize=0,
        linestyle="none",
    )
    finite = np.isfinite(ratio)
    if finite.any():
        span = np.nanpercentile(np.abs(ratio[finite] - 1.0), 95)
        span = max(0.15, min(1.5, span * 1.6))
        rax.set_ylim(1 - span, 1 + span)
    rax.set_ylabel("Data / MC", fontsize=9)
    rax.set_xlabel(label)
    rax.grid(axis="y", alpha=0.7)
    rax.set_axisbelow(True)

    fig.savefig(outpath, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    return outpath


def arm_overlay_plot(hists, outpath, sample_label):
    """Overlay the mass estimates on one selection.

    Confirms the three fit arms agree with each other and that the quality gate
    has not sculpted the peak. The raw pre-fit mass is drawn as a muted dashed
    reference: it is what the fits improve on, not a fourth arm.
    """
    all_gen = GEN_MC_CATEGORIES + [GEN_DATA]
    fig, ax = plt.subplots(figsize=(7.4, 4.6))

    raw_key = f"nominal_{RAW_MASS_VAR[0]}"
    if raw_key in hists:
        edges, vals, _ = project(hists[raw_key], all_gen)
        ax.step(
            edges,
            np.append(vals, vals[-1]),
            where="post",
            color=INK_MUTED,
            linewidth=1.6,
            linestyle="--",
            label=RAW_MASS_VAR[1],
        )

    for (var, name), color in zip(ARM_MASS_VARS, ARM_COLORS):
        key = f"nominal_{var}"
        if key not in hists:
            logger.warning("%s absent, skipping in overlay", key)
            continue
        edges, vals, _ = project(hists[key], all_gen)
        if vals.sum() <= 0:
            continue
        ax.step(
            edges,
            np.append(vals, vals[-1]),
            where="post",
            color=color,
            linewidth=2.0,
            label=name,
        )

    ax.set_xlabel(r"$m(\mu\mu K)$ [GeV]")
    ax.set_ylabel("Candidates / bin")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", alpha=0.7)
    ax.set_axisbelow(True)
    # Upper left is empty here: the B peak sits at 5.28, right of centre.
    ax.legend(loc="upper left", fontsize=9)
    ax.set_title(
        f"Mass estimates on one selection — {sample_label}",
        fontsize=10,
        color=INK,
        loc="left",
    )
    fig.savefig(outpath, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    return outpath


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--infile", required=True, help="Histmaker HDF5 output")
    parser.add_argument("-o", "--outfolder", default="plots", help="Output directory")
    parser.add_argument(
        "--variables",
        nargs="*",
        default=None,
        help="Variables to compare; default is a standard set",
    )
    args = parser.parse_args()

    style()
    os.makedirs(args.outfolder, exist_ok=True)
    samples = load(args.infile)
    data_hists = samples[DATA_SAMPLE]
    mc_hists = samples[MC_SAMPLE]

    variables = {v: v for v in args.variables} if args.variables else DEFAULT_VARIABLES
    written = []
    for var, label in variables.items():
        key = f"nominal_{var}"
        if key not in data_hists or key not in mc_hists:
            logger.warning("%s missing from one of the samples, skipping", key)
            continue
        path = os.path.join(args.outfolder, f"datamc_{var}.png")
        if data_mc_plot(data_hists[key], mc_hists[key], label, path):
            written.append(path)

    path = os.path.join(args.outfolder, "arm_overlay_data.png")
    written.append(arm_overlay_plot(data_hists, path, "data"))

    print(f"\nWrote {len(written)} plots to {args.outfolder}:")
    for path in written:
        print(f"  {path}")


if __name__ == "__main__":
    main()
