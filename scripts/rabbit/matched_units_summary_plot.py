#!/usr/bin/env python3
"""Summary figure for the A/e/M unit fix: old unit vs matched, per fit.

Four panels, one quantity each, fits on the y axis: sigma(POI) with the
"A, e, M fixed" floor marked, the impact of A+e+M on the POI, and the mean
posterior sigma of A_i and e_i in prior units. Colour is the arm (reference
without B+, old unit, matched); data fits are hatched. Run in the container.

    python scripts/rabbit/matched_units_summary_plot.py -o out.png --poi massShiftZ100MeV \
        --poiScale 100 --poiLabel mZ \
        --fit ref asimov "Z + J/psi (no B+)" fit1.hdf5 --fit old asimov "old unit" fit2.hdf5 ...
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, "rabbit")
sys.path.insert(0, "wums")
from rabbit import io_tools

ARM_COLOR = {"ref": "#8a8985", "old": "#eb6834", "matched": "#2a78d6"}
TEXT, TEXT2, GRID = "#0b0b0b", "#52514e", "#e6e6e3"

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--fit", nargs=4, action="append", required=True,
                metavar=("ARM", "KIND", "LABEL", "PATH"),
                help="ARM in ref/old/matched, KIND in asimov/data")
ap.add_argument("--poi", required=True)
ap.add_argument("--poiScale", type=float, required=True, help="MeV per unit of the POI")
ap.add_argument("--poiLabel", required=True)
ap.add_argument("--nEta", type=int, default=24)
ap.add_argument("--title", default="")
ap.add_argument("-o", "--output", required=True)
args = ap.parse_args()


def analyse(path):
    r = io_tools.get_fitresult(path)
    names = np.array([str(x) for x in r["parms"].get().axes["parms"]])
    cov = np.asarray(r["cov"].get().values())
    idx = {p: np.array([int(np.flatnonzero(names == f"{p}{i}")[0]) for i in range(args.nEta)])
           for p in "AeM"}
    sig = np.sqrt(np.diag(cov))
    out = dict(sigA=float(sig[idx["A"]].mean()), sige=float(sig[idx["e"]].mean()),
               sigma=np.nan, floor=np.nan, impact=np.nan)
    hit = np.flatnonzero(names == args.poi)
    if not hit.size:  # a reference fit without the POI (e.g. J/psi alone)
        return out
    ip = int(hit[0])
    fix = np.r_[idx["A"], idx["e"], idx["M"]]
    c = cov[ip, fix]
    cond = np.sqrt(cov[ip, ip] - c @ np.linalg.solve(cov[np.ix_(fix, fix)], c))
    s = sig[ip]
    out.update(sigma=args.poiScale * s, floor=args.poiScale * cond,
               impact=args.poiScale * np.sqrt(max(s * s - cond * cond, 0.0)))
    return out


fits = [(arm, kind, label, analyse(path)) for arm, kind, label, path in args.fit]

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

panels = [
    ("sigma", f"σ({args.poiLabel}) [MeV]", "{:.4f}" if args.poiScale < 10 else "{:.3f}"),
    ("impact", f"impact of A+e+M on {args.poiLabel} [MeV]", "{:.4f}" if args.poiScale < 10 else "{:.3f}"),
    ("sigA", "mean σ(A_i) [prior units]", "{:.3f}"),
    ("sige", "mean σ(e_i) [prior units]", "{:.3f}"),
]
n = len(fits)
fig, axs = plt.subplots(1, 4, figsize=(15, 0.55 * n + 1.9), sharey=True)
y = np.arange(n)[::-1]
for ax, (key, title, fmt) in zip(axs, panels):
    vals = np.array([f[3][key] for f in fits])
    for yi, (arm, kind, label, o), v in zip(y, fits, vals):
        if not np.isfinite(v):
            ax.text(0, yi, " — (no POI)", va="center", ha="left", fontsize=8.5, color=TEXT2)
            continue
        ax.barh(yi, v, height=0.62, color=ARM_COLOR[arm], edgecolor="#fcfcfb", linewidth=2,
                hatch="//" if kind == "data" else None)
        ax.text(v, yi, " " + fmt.format(v), va="center", ha="left", fontsize=8.5, color=TEXT)
        if key == "sigma":
            ax.plot([o["floor"]] * 2, [yi - 0.36, yi + 0.36], color=TEXT, lw=2)
    ax.set_xlim(0, np.nanmax(vals) * 1.35)
    ax.set_title(title, fontsize=10.5, color=TEXT)
    ax.grid(axis="x", color=GRID)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=TEXT2, labelsize=8.5)
axs[0].set_yticks(y)
axs[0].set_yticklabels([f"{label} ({kind})" for _, kind, label, _ in fits], fontsize=9, color=TEXT)

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

handles = [Patch(color=ARM_COLOR[a], label=l) for a, l in
           (("ref", "no B+"), ("old", "B+, old unit (massfit widths)"), ("matched", "B+, matched unit"))
           if any(f[0] == a for f in fits)]
handles += [Patch(facecolor="white", edgecolor=TEXT2, hatch="//", label="DATA (mB blinded)")] \
    if any(f[1] == "data" for f in fits) else []
handles += [Line2D([], [], color=TEXT, lw=2, label=f"σ({args.poiLabel}) with A, e, M fixed")]
fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, fontsize=9,
           bbox_to_anchor=(0.5, -0.01))
if args.title:
    fig.suptitle(args.title, fontsize=12, color=TEXT)
fig.tight_layout(rect=(0, 0.08, 1, 0.97 if args.title else 1))
fig.savefig(args.output, dpi=130)
fig.savefig(args.output.rsplit(".", 1)[0] + ".pdf")
print("wrote", args.output)
