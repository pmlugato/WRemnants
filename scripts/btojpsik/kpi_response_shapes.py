"""Response-reweight shapes for the muon and kaon/pion ONNX models.

Pure network probes: load each ONNX with onnxruntime and evaluate
``log W`` on a controlled grid. No data, no histmaker, no fit -- which is the
point, because the failure this is looking for (a model asked to extrapolate
outside its training range) makes a template variation LARGER, and a larger
variation always reads as a SMALLER uncertainty in prior units. An
out-of-domain model therefore looks like improved sensitivity in any fit
output, and cannot be caught by inspecting one.

Produces, for openspec/changes/add-kpi-onnx-reweights section 2:
  2.1  slope d log W / d(delta r_kappa) against pt, per model and class
  2.2  the same against eta, and against the response residual y_raw, to check
       the single-slice probes are representative rather than a lucky point
  2.3  full log W against delta r_kappa at low and high pt -- the shape, not
       just the slope, since a fit sees the whole variation and non-linearity
       is what the conditioning clamp was hiding

Example:
    python scripts/btojpsik/kpi_response_shapes.py -o ~/public_html/mz/kpi_260908
"""

import argparse
import os
import sys

import matplotlib as mpl
import numpy as np
import onnxruntime as ort

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plotting"))
from btojpsik_cvhnano_datamc import INK, INK_MUTED, SURFACE, style  # noqa: E402

DATA_DIR = "wremnants-data/data/calibration"
STEM = "shift_smear_reweight_mlp_factored_combined"

# (label, file suffix, class code, colour, linestyle). Codes are the alphabet
# of the model that reads them -- Muon_genPartFlav for the muon family, PDG
# ids for kpi -- and they are NOT interchangeable: 321 through the muon
# family's internal remap becomes 1, its prompt class.
SERIES = [
    ("muon model, J/psi leg (443)", "_nosmearing", 443, "#2a78d6", "-"),
    ("muon model, prompt (1)", "_nosmearing", 1, "#2a78d6", "--"),
    ("muon model, tau (15)", "_nosmearing", 15, "#2a78d6", ":"),
    ("kpi model, kaon (321)", "_kpi", 321, "#eb6834", "-"),
    ("kpi model, pion (211)", "_kpi", 211, "#eb6834", "--"),
    ("kpi model, muon (13)", "_kpi", 13, "#eb6834", ":"),
]

# Recorded usable windows, from design.md. Shaded on the pt plot.
WINDOWS = {"_nosmearing": (2.0, None), "_kpi": (1.0, 15.0)}

_SESS = {}


def _logx_ticks(ax, vals):
    """Readable ticks on a log x axis.

    matplotlib's default LogLocator prints every minor decade multiple, which
    on a range of one or two decades collides into an unreadable smear. Set
    them explicitly instead.
    """
    from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter
    ax.xaxis.set_major_locator(FixedLocator(vals))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda v, _: f"{v:g}"))


def sess(suffix):
    if suffix not in _SESS:
        _SESS[suffix] = ort.InferenceSession(
            f"{DATA_DIR}/{STEM}{suffix}.onnx", providers=["CPUExecutionProvider"]
        )
    return _SESS[suffix]


def logw(suffix, pt, eta, code, u, y=(0.0, 0.0, 0.0), phi=0.3, q=1.0):
    """log W for an array of delta r_kappa values at one kinematic point."""
    u = np.atleast_1d(np.asarray(u, dtype=np.float32))
    c = np.array(
        [[np.log(pt), q, np.arctan(np.sinh(eta)), np.sin(phi), np.cos(phi), float(code)]],
        dtype=np.float32,
    )
    y_raw = np.asarray([y], dtype=np.float32)
    U = np.zeros((1, u.size, 3), dtype=np.float32)
    U[0, :, 0] = u
    S = np.zeros((1, u.size, 3), dtype=np.float32)
    return sess(suffix).run(
        ["d"], {"y_raw": y_raw, "c_raw": c, "u": U, "sigma": S}
    )[0][0]


def slope(suffix, pt, eta, code, h=1e-5, **kw):
    d = logw(suffix, pt, eta, code, np.array([h, -h], dtype=np.float32), **kw)
    return (d[0] - d[1]) / (2 * h)


def _finish(fig, ax, out, xlabel, ylabel, title):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote", out)


def fig_slope_vs_pt(outdir):
    pts = np.geomspace(0.3, 60.0, 90)
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    for label, suf, code, col, ls in SERIES:
        ax.plot(pts, [slope(suf, p, 0.6, code) for p in pts],
                color=col, linestyle=ls, linewidth=1.6, label=label)
    ax.axhline(0.0, color=INK_MUTED, linewidth=0.9)
    for suf, col in (("_nosmearing", "#2a78d6"), ("_kpi", "#eb6834")):
        lo, hi = WINDOWS[suf]
        ax.axvspan(lo, hi if hi else pts[-1], color=col, alpha=0.05)
    ax.set_xscale("log")
    ax.set_xlim(pts[0], pts[-1])
    _logx_ticks(ax, [0.3, 0.5, 1, 2, 5, 10, 20, 50])
    ax.legend(fontsize=8, ncol=2, frameon=False)
    _finish(fig, ax, os.path.join(outdir, "slope_vs_pt.png"),
            r"gen $p_T$ [GeV]", r"$\mathrm{d}\log W / \mathrm{d}(\delta r_\kappa)$",
            r"Response slope at $\eta$ = 0.6, $q$ = +1, $y_{raw}$ = 0"
            "   (shaded: each model's recorded window)")


def fig_slope_vs_eta(outdir):
    etas = np.linspace(-2.3, 2.3, 47)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), sharey=True)
    for ax, pt in zip(axes, (2.0, 10.0)):
        for label, suf, code, col, ls in SERIES:
            ax.plot(etas, [slope(suf, pt, e, code) for e in etas],
                    color=col, linestyle=ls, linewidth=1.5, label=label)
        ax.axhline(0.0, color=INK_MUTED, linewidth=0.9)
        ax.set_title(rf"gen $p_T$ = {pt:g} GeV", fontsize=10, color=INK, loc="left")
        ax.set_xlabel(r"gen $\eta$")
        ax.grid(alpha=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel(r"$\mathrm{d}\log W / \mathrm{d}(\delta r_\kappa)$")
    axes[0].legend(fontsize=7.5, frameon=False)
    out = os.path.join(outdir, "slope_vs_eta.png")
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote", out)


def fig_slope_vs_yraw(outdir):
    """Every other probe here sits at y_raw = 0. This is the check that the
    conclusions are not an artefact of that one slice."""
    # r_kappa residual in units of each model's own target_std, so the two are
    # compared over the same part of their respective response distributions.
    zs = np.linspace(-3.0, 3.0, 41)
    std = {"_nosmearing": 0.016529, "_kpi": 0.012567}
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    for label, suf, code, col, ls in SERIES:
        ax.plot(zs, [slope(suf, 2.0, 0.6, code, y=(z * std[suf], 0.0, 0.0))
                     for z in zs],
                color=col, linestyle=ls, linewidth=1.5, label=label)
    ax.axhline(0.0, color=INK_MUTED, linewidth=0.9)
    ax.axvline(0.0, color=INK_MUTED, linewidth=0.9, linestyle=":")
    ax.legend(fontsize=8, ncol=2, frameon=False)
    _finish(fig, ax, os.path.join(outdir, "slope_vs_yraw.png"),
            r"response residual $r_\kappa$ / target_std",
            r"$\mathrm{d}\log W / \mathrm{d}(\delta r_\kappa)$",
            r"Slope against where the muon sits in its own response, "
            r"$p_T$ = 2 GeV, $\eta$ = 0.6")


def fig_cliff_vs_yraw(outdir):
    """The decisive one: slope against pt at several response residuals, for
    the smeared-trained model against the un-smeared one.

    Every other probe here sits at y_raw = 0, and fig_slope_vs_yraw shows that
    is the WORST place to sit -- the slope is a steep odd-ish function of y_raw
    with amplitude ~250 and a zero crossing right there. So a pt trend measured
    at y_raw = 0 is measuring the residual left over after the large offset
    cancels, not the response.
    """
    pts = np.geomspace(0.7, 8.0, 50)
    zs = [-2.0, -1.0, 0.0, 1.0, 2.0]
    cols = plt.cm.coolwarm(np.linspace(0.05, 0.95, len(zs)))
    std = 0.016529
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2), sharey=True)
    for ax, suf, name in (
        (axes[0], "_smearing", "muon model, SMEARED-trained"),
        (axes[1], "_nosmearing", "muon model, un-smeared (correct here)"),
    ):
        for z, c in zip(zs, cols):
            ax.plot(pts,
                    [slope(suf, p, 0.6, 443, y=(z * std, 0.0, 0.0)) for p in pts],
                    color=c, linewidth=1.7,
                    label=rf"$r_\kappa$/std = {z:+.0f}")
        ax.axhline(0.0, color=INK_MUTED, linewidth=0.9)
        ax.axvline(1.8, color="#b4232a", linewidth=1.2, linestyle="--")
        ax.set_xscale("log")
        ax.set_xlim(pts[0], pts[-1])
        _logx_ticks(ax, [0.7, 1, 1.5, 2, 3, 5, 8])
        ax.set_xlabel(r"gen $p_T$ [GeV]")
        ax.set_title(name, fontsize=10, color=INK, loc="left")
        ax.grid(alpha=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel(r"$\mathrm{d}\log W / \mathrm{d}(\delta r_\kappa)$")
    axes[0].legend(fontsize=8, frameon=False, ncol=2)
    out = os.path.join(outdir, "cliff_vs_yraw.png")
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote", out)


def fig_logw_shape(outdir):
    us = np.linspace(-1.5e-3, 1.5e-3, 121).astype(np.float32)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    for ax, pt in zip(axes, (1.5, 10.0)):
        for label, suf, code, col, ls in SERIES:
            ax.plot(us * 1e3, logw(suf, pt, 0.6, code, us),
                    color=col, linestyle=ls, linewidth=1.5, label=label)
        ax.axhline(0.0, color=INK_MUTED, linewidth=0.9)
        ax.axvline(0.0, color=INK_MUTED, linewidth=0.9, linestyle=":")
        ax.set_title(rf"gen $p_T$ = {pt:g} GeV", fontsize=10, color=INK, loc="left")
        ax.set_xlabel(r"$\delta r_\kappa$  [$\times 10^{-3}$]")
        ax.grid(alpha=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel(r"$\log W$")
    axes[0].legend(fontsize=7.5, frameon=False)
    out = os.path.join(outdir, "logw_shape.png")
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote", out)


def report():
    """The numbers behind the figures, for the record."""
    print("\nSanity: log W at delta r_kappa = 0 must be exactly 0")
    bad = 0
    for _, suf, code, _, _ in SERIES:
        for pt in (0.5, 2.0, 10.0, 45.0):
            v = logw(suf, pt, 0.6, code, np.array([0.0], dtype=np.float32))[0]
            bad += v != 0.0
    print(f"  {'PASS' if not bad else f'FAIL ({bad})'}")

    print("\nSlope table, eta = 0.6, q = +1, y_raw = 0")
    hdr = f"{'pT':>6} | " + " ".join(f"{l.split(', ')[1]:>12}" for l, *_ in SERIES)
    print(hdr)
    for pt in (0.3, 0.5, 1.0, 1.5, 1.95, 2.5, 4.0, 5.0, 10.0, 15.0, 25.0, 45.0):
        row = f"{pt:6.2f} | " + " ".join(
            f"{slope(suf, pt, 0.6, code):12.1f}" for _, suf, code, _, _ in SERIES)
        print(row)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-o", "--outfolder", required=True)
    args = p.parse_args()
    style()
    os.makedirs(args.outfolder, exist_ok=True)
    fig_slope_vs_pt(args.outfolder)
    fig_slope_vs_eta(args.outfolder)
    fig_slope_vs_yraw(args.outfolder)
    fig_cliff_vs_yraw(args.outfolder)
    fig_logw_shape(args.outfolder)
    report()


if __name__ == "__main__":
    main()
