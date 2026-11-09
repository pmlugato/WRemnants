"""Cutflow waterfall for the CVH NanoAOD B+ histmaker.

Parses the text table `--saveCutflow` writes and draws the surviving fraction
per selection, data and simulation side by side. Two panels because the two
tables answer different questions: an event with several candidates survives a
cut that kills most of them, so the event table understates every per-candidate
cut and the candidate table is the one that says what a cut costs.

Example:

    python scripts/btojpsik/cvhnano_cutflow_fig.py \
        --cutflow .../cutflow_cutvars_nom.txt -o ~/public_html/mz/cutvars/
"""

import argparse
import os
import re
import sys

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plotting"))
from btojpsik_cvhnano_datamc import COLOR_SIGNAL, INK, INK_MUTED, SURFACE, style  # noqa: E402

COLOR_MC = "#eb6834"
ROW = re.compile(r"^(.{1,42}?)\s{2,}([\d.eE+-]+)\s+\(([\d.]+)\)\s+([\d.eE+-]+)\s+\(([\d.]+)\)\s*$")


def parse(path):
    """Return [(title, [(name, frac_data, frac_mc), ...]), ...] per table."""
    tables, cur, title = [], [], None
    for line in open(path, encoding="utf-8"):
        if line.startswith("Cutflow,"):
            if cur:
                tables.append((title, cur))
            title, cur = line.strip(), []
            continue
        m = ROW.match(line.rstrip())
        if m:
            cur.append((m.group(1).strip(), float(m.group(3)), float(m.group(5))))
    if cur:
        tables.append((title, cur))
    return tables


def panel(ax, rows, label):
    names = [r[0] for r in rows]
    y = np.arange(len(names))
    ax.barh(y - 0.2, [r[1] for r in rows], height=0.38, color=COLOR_SIGNAL,
            label="Data")
    ax.barh(y + 0.2, [r[2] for r in rows], height=0.38, color=COLOR_MC,
            label="Simulation")
    # Labels go past the END of the longer bar. Putting them before the shorter
    # one collides with the tick labels whenever a cut bites hard, which is
    # exactly the row worth reading.
    for i, r in enumerate(rows):
        ax.text(max(r[1], r[2]) + 0.012, i, f"{r[1]:.3f} / {r[2]:.3f}",
                ha="left", va="center", fontsize=7.5, color=INK_MUTED)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.20)
    ax.set_xlabel(label, fontsize=9)
    ax.grid(axis="x", alpha=0.7)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", fontsize=8)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cutflow", required=True)
    p.add_argument("-o", "--outfolder", required=True)
    p.add_argument("--tag", default="nominal")
    args = p.parse_args()

    style()
    os.makedirs(args.outfolder, exist_ok=True)
    tables = parse(args.cutflow)
    if not tables:
        raise SystemExit(f"no cutflow rows parsed from {args.cutflow}")

    n = len(tables)
    fig, axes = plt.subplots(1, n, figsize=(6.6 * n, 0.30 * max(len(t[1]) for t in tables) + 2.2))
    axes = np.atleast_1d(axes)
    for ax, (title, rows) in zip(axes, tables):
        panel(ax, rows, "events surviving / total" if "EVENTS" in title
              else "candidates surviving / total")
        ax.set_title("Events with >= 1 candidate" if "EVENTS" in title
                     else "Candidates", fontsize=10, color=INK, loc="left")
    fig.suptitle(f"B$^+\\to$J/$\\psi$K$^+$ CVH NanoAOD cutflow — {args.tag}",
                 fontsize=11, color=INK, x=0.005, ha="left")
    out = os.path.join(args.outfolder, f"cutflow_{args.tag}.png")
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
