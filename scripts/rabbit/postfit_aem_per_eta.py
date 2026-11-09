#!/usr/bin/env python3
"""Postfit A_i / e_i / M_i (value +- sigma, prior units) per eta group, several fits overlaid.

    python scripts/rabbit/postfit_aem_per_eta.py -o out.png --fit LABEL COLOR fit.hdf5 [--fit ...]
"""
import argparse, sys
import numpy as np
sys.path.insert(0, "rabbit"); sys.path.insert(0, "wums")
from rabbit import io_tools

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--fit", nargs=3, action="append", required=True, metavar=("LABEL", "COLOR", "PATH"))
ap.add_argument("--nEta", type=int, default=24)
ap.add_argument("--units", nargs=3, type=float, default=[5e-3, 5e-2, 1e-4], help="one unit of A, e, M")
ap.add_argument("--title", default="")
ap.add_argument("-o", "--output", required=True)
a = ap.parse_args()

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
TEXT, TEXT2, GRID = "#0b0b0b", "#52514e", "#e6e6e3"
eta = np.linspace(-2.4, 2.4, a.nEta + 1); etac = 0.5 * (eta[1:] + eta[:-1])
fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
n = len(a.fit); off = (np.arange(n) - (n - 1) / 2) * 0.045
for k, (label, color, path) in enumerate(a.fit):
    h = io_tools.get_fitresult(path)["parms"].get()
    names = np.array([str(x) for x in h.axes["parms"]]); v = h.values(); s = np.sqrt(h.variances())
    for ax, p in zip(axs, "AeM"):
        ii = [int(np.flatnonzero(names == f"{p}{i}")[0]) for i in range(a.nEta)]
        ax.errorbar(etac + off[k], v[ii], s[ii], fmt="o", ms=4.5, lw=1.4, color=color, label=label)
for ax, p, u in zip(axs, "AeM", a.units):
    ax.axhline(0, color=TEXT2, lw=0.8); ax.grid(color=GRID); ax.set_axisbelow(True)
    ax.set_ylabel(f"postfit {p}_i [units of {u:g}]", color=TEXT)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axs[0].legend(frameon=False, fontsize=9, ncol=n, loc="upper center")
axs[-1].set_xlabel("eta group centre", color=TEXT)
if a.title: fig.suptitle(a.title, fontsize=11.5, color=TEXT)
fig.tight_layout(); fig.savefig(a.output, dpi=130); fig.savefig(a.output.rsplit(".", 1)[0] + ".pdf"); print("wrote", a.output)
