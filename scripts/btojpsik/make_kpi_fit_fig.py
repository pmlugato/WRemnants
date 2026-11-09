"""Section 5 and 6 results as a figure, for the deck.

Left: the standalone B+ Asimov, sigma in prior units per configuration.
Right: the joint fit, sigma(mZ) relative to the Z + J/psi reference.

Both are plotted as ratios rather than absolutes because the absolute values
depend on the cell binning and the reference chosen, while the comparison
between configurations -- the thing the study is for -- does not.
"""
import os
import sys

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plotting"))
from btojpsik_cvhnano_datamc import INK, INK_MUTED, SURFACE, style  # noqa: E402

CFG = ["muon_rec", "muon_cor", "kpi", "perleg"]
NICE = ["muon\n(recorded)", "muon\n(smearing fixed)", "kpi", "kpi + per-leg\n(config D)"]
COLS = ["#8a8a88", "#2a78d6", "#eb6834", "#4a3aa7"]

# Standalone Asimov, mean sigma in prior units (openspec 5.3)
SOLO = {"A": [0.9531, 0.9476, 0.9570, 0.9560],
        "e": [0.6845, 0.6720, 0.7191, 0.7188],
        "rho": [0.2289, 0.2419, 0.1995, 0.2008]}
# Joint Asimov, sigma(mZ) ratio to the Z + J/psi reference (openspec 6.2)
JOINT = [0.9968, 0.9965, 0.9975, 0.9975]


def main(outdir):
    style()
    os.makedirs(outdir, exist_ok=True)
    x = np.arange(len(CFG))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.2))

    w = 0.36
    a1.bar(x - w / 2, SOLO["e"], w, color=COLS, edgecolor=SURFACE,
           linewidth=1.0, label=r"$\sigma(e)$")
    a1.bar(x + w / 2, SOLO["rho"], w, color=COLS, alpha=0.45,
           edgecolor=SURFACE, linewidth=1.0, label=r"$\rho(A,e)$")
    for xi, (e, r) in enumerate(zip(SOLO["e"], SOLO["rho"])):
        a1.text(xi - w / 2, e + 0.012, f"{e:.3f}", ha="center", fontsize=8.5,
                color=INK)
        a1.text(xi + w / 2, r + 0.012, f"{r:.3f}", ha="center", fontsize=8.5,
                color=INK_MUTED)
    a1.set_xticks(x)
    a1.set_xticklabels(NICE, fontsize=8.5)
    a1.set_ylim(0, 0.85)
    a1.set_ylabel("prior units")
    a1.set_title("Standalone B$^{+}$ Asimov: kpi is slightly weaker on e,\n"
                 "slightly better on the A-e degeneracy",
                 fontsize=9.5, color=INK, loc="left")
    a1.legend(fontsize=8.5, frameon=False, loc="upper right")
    a1.grid(axis="y", alpha=0.6)
    a1.set_axisbelow(True)

    gain = [(1.0 - g) * 100 for g in JOINT]
    a2.bar(x, gain, 0.5, color=COLS, edgecolor=SURFACE, linewidth=1.0)
    for xi, g in enumerate(gain):
        a2.text(xi, g + 0.006, f"{g:.2f}%", ha="center", fontsize=9, color=INK)
    a2.set_xticks(x)
    a2.set_xticklabels(NICE, fontsize=8.5)
    a2.set_ylim(0, 0.46)
    a2.set_ylabel(r"improvement in $\sigma(m_Z)$ [%]")
    a2.set_title("Joint Z + J/$\\psi$ + B$^{+}$ Asimov: the 0.3% gain\n"
                 "survives, and kpi shaves it slightly",
                 fontsize=9.5, color=INK, loc="left")
    a2.grid(axis="y", alpha=0.6)
    a2.set_axisbelow(True)

    out = os.path.join(outdir, "kpi_fit_summary.png")
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[sys.argv.index("-o") + 1])
