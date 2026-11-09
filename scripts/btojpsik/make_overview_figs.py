#!/usr/bin/env python3
"""Three figures for the B+ -> J/psi K overview deck.

  A  model ladder: chi2/ndf against free parameters, physics vs flexibility
  B  what more signal simulation buys, against the data-statistical floor
  C  the momentum-scale lever arm each channel of the combined fit covers
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = os.path.expanduser("~/public_html/slides/assets_btojpsik_overview")
os.makedirs(OUT, exist_ok=True)
INK, GREY, ACC = "#1F1F1F", "#7F7F7F", "#C0392B"


# ------------------------------------------------------------------ A
def model_ladder():
    flex = [(4, 531.75, "poly1"), (5, 53.59, "poly2"), (6, 51.01, "poly3"),
            (7, 48.43, "poly4"), (9, 22.40, "poly6")]
    phys = [(5, 53.59, ""), (6, 44.16, "+asym res"), (8, 26.43, "+3 templates"),
            (9, 10.10, "both"), (10, 10.21, ""), (11, 8.31, "")]
    f, ax = plt.subplots(figsize=(9.8, 4.9))
    x, y, _ = zip(*flex)
    ax.plot(x, y, "o-", color=GREY, ms=7, lw=1.8,
            label="more analytic flexibility (exp·poly$N$ + Gaussian resolution)")
    x2, y2, _ = zip(*phys)
    ax.plot(x2, y2, "s-", color=ACC, ms=7, lw=2.0,
            label="physics ingredients (asymmetric resolution, $B^0$/$J\\psi\\pi$/$J\\psi KX$)")
    for xx, yy, lab in flex:
        if lab in ("poly2", "poly6"):
            ax.annotate(lab, (xx, yy), textcoords="offset points", xytext=(4, 8),
                        fontsize=9, color=GREY)
    for xx, yy, lab in phys:
        if lab:
            ax.annotate(lab, (xx, yy), textcoords="offset points", xytext=(-52, -18),
                        fontsize=9, color=ACC)
    ax.annotate("", xy=(9, 10.10), xytext=(9, 22.40),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.3))
    ax.text(9.15, 15.0, "same 9 parameters,\n2.2× better", fontsize=9.5, color=INK)
    ax.set_yscale("log")
    ax.set_xlabel("free parameters in the model")
    ax.set_ylabel("$\\chi^2$/ndf on the data mass spectrum")
    ax.set_xlim(3.4, 11.8)
    ax.legend(frameon=False, fontsize=9.5, loc="upper right")
    ax.grid(alpha=0.15)
    p = os.path.join(OUT, "model_ladder.png")
    f.savefig(p, dpi=150, bbox_inches="tight"); plt.close(f)
    print("MARKER wrote", p)


# ------------------------------------------------------------------ B
def mc_scaling():
    floor, mc = 0.4280, 0.5860          # measured on the converged Asimov pair
    N = np.logspace(0, 2.2, 400)
    s = np.sqrt(floor ** 2 + mc ** 2 / N)
    s0 = np.sqrt(floor ** 2 + mc ** 2)
    f, ax = plt.subplots(figsize=(9.8, 4.7))
    ax.plot(N, 100 * s / s0, "-", color=ACC, lw=2.2)
    ax.axhline(100 * floor / s0, color=INK, ls="--", lw=1.3)
    ax.text(1.05, 100 * floor / s0 - 3.4,
            "data-statistical floor — no amount of simulation goes past it",
            fontsize=9.5, color=INK)
    for n, lab in ((1, "today — 0.90 fb$^{-1}$ of\nsignal simulation"), (6, "6×\n4.3M excl."),
                   (18.7, "data parity\n13.4M excl.")):
        v = 100 * np.sqrt(floor ** 2 + mc ** 2 / n) / s0
        ax.plot([n], [v], "o", color=INK, ms=7)
        ax.annotate(lab, (n, v), textcoords="offset points", xytext=(6, 10),
                    fontsize=9, color=INK)
    ax.set_xscale("log")
    ax.set_xlabel("signal simulation, relative to today")
    ax.set_ylabel("A/e/M uncertainty  [% of today]")
    ax.set_xticks([1, 2, 4, 6, 10, 20, 50, 100])
    ax.set_xticklabels(["1", "2", "4", "6", "10", "20", "50", "100"])
    ax.set_ylim(52, 103)
    ax.grid(alpha=0.15)
    p = os.path.join(OUT, "mc_scaling.png")
    f.savefig(p, dpi=150, bbox_inches="tight"); plt.close(f)
    print("MARKER wrote", p)


# ------------------------------------------------------------------ C
def lever_arm():
    """dk/k = A - e*k + q*M/k : each channel covers a different range in k = 1/pT."""
    k = np.logspace(np.log10(0.012), np.log10(1.1), 500)
    f, (a1, a2) = plt.subplots(2, 1, figsize=(9.8, 5.6), sharex=True,
                               gridspec_kw=dict(height_ratios=[1.9, 1], hspace=0.08))
    a1.plot(k, 1e4 * (1e-3 * np.ones_like(k)), color=INK, lw=2.0, label="$A$ only  ($10^{-3}$)")
    a1.plot(k, 1e4 * (-1e-2 * k), color=ACC, lw=2.0, label="$e$ only  ($10^{-2}$)")
    a1.plot(k, 1e4 * (1e-4 / k), color="#2E86C1", lw=2.0, ls="--", label="$M$ only  ($10^{-4}$)")
    a1.axhline(0, color=GREY, lw=0.7)
    a1.set_ylabel("$\\delta k/k$   [$10^{-4}$]")
    a1.set_ylim(-26, 26)
    a1.legend(frameon=False, fontsize=9.5, ncol=3, loc="upper left")

    chans = [
        ("Z $\\to\\mu\\mu$",        25.0, 60.0, "#2E86C1"),
        ("J/$\\psi\\to\\mu\\mu$",    3.0, 20.0, "#27AE60"),
        ("D$^0\\to$K$\\pi$",         1.0, 10.0, "#8E44AD"),
        ("B$^+\\to$J/$\\psi$K, kaon", 1.0, 10.0, ACC),
    ]
    for i, (lab, lo, hi, col) in enumerate(chans):
        a2.plot([1 / hi, 1 / lo], [i, i], lw=7, color=col, solid_capstyle="butt")
        a2.text(1 / lo * 1.10, i, lab, va="center", fontsize=10, color=INK)
    a2.set_ylim(-0.7, len(chans) - 0.1)
    a2.set_yticks([])
    a2.set_xlabel("curvature  $k = 1/p_T$   [GeV$^{-1}$]")
    for a in (a1, a2):
        a.set_xscale("log")
        a.set_xlim(0.012, 1.9)
        a.grid(alpha=0.15, axis="x")
    p = os.path.join(OUT, "lever_arm.png")
    f.savefig(p, dpi=150, bbox_inches="tight"); plt.close(f)
    print("MARKER wrote", p)


model_ladder()
mc_scaling()
lever_arm()
