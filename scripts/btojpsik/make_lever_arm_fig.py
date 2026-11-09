#!/usr/bin/env python3
"""The A-flatness test: why the unclamped sensitivity is not usable.

Left  -- the induced peak shift per bachelor-pT band. A is a pT-INDEPENDENT
         scale, so its curve must be flat. It is, clamped; it is destroyed when
         the network extrapolates below its 2 GeV training floor.
Right -- what the three Asimov fits report for sigma(e) and rho(A,e).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = os.path.expanduser("~/public_html/slides/assets_btojpsik_overview")
INK, ACC, BLU, GREY = "#1F1F1F", "#C0392B", "#2E86C1", "#7F7F7F"

bands = [(0.1, 0.3), (0.3, 0.5), (0.5, 0.8), (0.8, 1.0), (1.0, 1.5),
         (1.5, 2.0), (2.0, 3.0), (3.0, 5.0), (5.0, 10.0)]
Ac = [-0.039, -0.038, -0.039, -0.042, -0.052, -0.057, -0.057, -0.057, -0.056]
An = [-0.380, 0.003, -0.005, 0.012, 0.016, -0.034, -0.057, -0.057, -0.056]
ec = [1.162, 0.767, 0.552, 0.412, 0.344, 0.270, 0.194, 0.121, 0.070]
en = [1.544, 0.021, 0.326, -0.050, -0.117, 0.155, 0.194, 0.121, 0.070]
x = np.array([np.sqrt(a * b) for a, b in bands])

f, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 4.4))
a1.plot(x, Ac, "o-", color=INK, ms=6, lw=2.0, label="$A$, clamped at 2 GeV")
a1.plot(x, An, "o--", color=INK, ms=5, lw=1.4, alpha=0.5, mfc="none",
        label="$A$, network extrapolated")
a1.plot(x, ec, "s-", color=ACC, ms=6, lw=2.0, label="$e$, clamped")
a1.plot(x, en, "s--", color=ACC, ms=5, lw=1.4, alpha=0.5, mfc="none",
        label="$e$, extrapolated")
a1.axvline(2.0, color=GREY, ls=":", lw=1.4)
a1.text(2.15, -0.30, "training\nfloor", fontsize=9, color=GREY)
a1.axhline(0, color=GREY, lw=0.7)
a1.set_xscale("log")
a1.set_xlabel("bachelor $p_T$  [GeV]")
a1.set_ylabel("induced peak shift  [MeV]")
a1.set_title("$A$ must be flat in $p_T$ — and is not, unclamped", fontsize=10)
a1.legend(frameon=False, fontsize=8.5, ncol=2, loc="upper center")
a1.grid(alpha=0.15)

labs = ["$p_T>1.0$\nclamped", "$p_T>0.1$\nclamped", "$p_T>0.1$\nclamp OFF"]
se = [0.343, 0.285, 0.074]
rho = [0.586, 0.540, 0.908]
i = np.arange(3)
a2.bar(i - 0.19, se, 0.36, color=ACC, label=r"$\sigma(e)$  [prior units]")
a2.bar(i + 0.19, rho, 0.36, color=BLU, label=r"$\rho(A,e)$")
for k in i:
    a2.text(k - 0.19, se[k] + 0.02, f"{se[k]:.3f}", ha="center", fontsize=9, color=ACC)
    a2.text(k + 0.19, rho[k] + 0.02, f"{rho[k]:.3f}", ha="center", fontsize=9, color=BLU)
a2.set_xticks(i); a2.set_xticklabels(labs, fontsize=9)
a2.set_ylim(0, 1.08)
a2.legend(frameon=False, fontsize=9)
a2.set_title("The two treatments disagree 4.6× on $\\sigma(e)$", fontsize=10)
a2.grid(alpha=0.15, axis="y")
p = os.path.join(OUT, "lever_arm_clamp_test.png")
f.savefig(p, dpi=150, bbox_inches="tight")
print("MARKER wrote", p)
