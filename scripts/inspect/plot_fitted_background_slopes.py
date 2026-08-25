#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from rabbit import io_tools


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Check whether fitted exponential background slopes can be shared across "
            "pt bins for each (eta1, eta2) cell."
        )
    )
    parser.add_argument("fitresult", help="Rabbit fitresults HDF5 file")
    parser.add_argument("-o", "--outdir", default="fitted_background_slopes")
    parser.add_argument("--process", default="background_JPsi")
    parser.add_argument("--eta1", type=int, help="eta1 bin for the detailed plot")
    parser.add_argument("--eta2", type=int, help="eta2 bin for the detailed plot")
    parser.add_argument(
        "--useCovariance",
        action="store_true",
        help="Use the full fitted covariance for compatibility tests (slower and heavier)",
    )
    parser.add_argument("--top", type=int, default=12)
    return parser.parse_args()


def load_slopes(path, process, use_covariance):
    result = io_tools.get_fitresult(path)
    hparams = result["parms"].get()
    names = np.asarray(hparams.axes["parms"], dtype=str)
    values = np.asarray(hparams.values())
    variances = np.asarray(hparams.variances())
    pattern = re.compile(
        rf"^slope_{re.escape(process)}_"
        r"eta1(?P<eta1>\d+)_eta2(?P<eta2>\d+)_"
        r"pt1(?P<pt1>\d+)_pt2(?P<pt2>\d+)$"
    )

    slopes = []
    for index, (name, value, variance) in enumerate(zip(names, values, variances)):
        match = pattern.match(name)
        if not match:
            continue
        slopes.append(
            {
                "index": index,
                "name": name,
                "value": float(value),
                "variance": float(variance),
                **{key: int(value) for key, value in match.groupdict().items()},
            }
        )

    if not slopes:
        raise RuntimeError(
            f"No slope parameters found for process '{process}'. "
            "Check --process against the fitted parameter names."
        )

    covariance = None
    if use_covariance:
        covariance = np.asarray(result["cov"].get().values())

    return slopes, covariance


def compatibility(group, covariance):
    values = np.asarray([entry["value"] for entry in group])
    indices = np.asarray([entry["index"] for entry in group])
    variances = np.asarray([entry["variance"] for entry in group])
    finite = np.isfinite(values) & np.isfinite(variances) & (variances > 0.0)
    values = values[finite]
    indices = indices[finite]
    variances = variances[finite]

    if len(values) < 2:
        return np.nan, np.nan, np.nan, len(values)

    if covariance is None:
        inverse = np.diag(1.0 / variances)
    else:
        subcov = covariance[np.ix_(indices, indices)]
        inverse = np.linalg.pinv(subcov, hermitian=True)

    ones = np.ones(len(values))
    denominator = ones @ inverse @ ones
    if not np.isfinite(denominator) or denominator <= 0.0:
        return np.nan, np.nan, float(np.std(values)), len(values)

    mean = float((ones @ inverse @ values) / denominator)
    residual = values - mean
    chi2 = float(residual @ inverse @ residual)
    return mean, chi2 / (len(values) - 1), float(np.std(values)), len(values)


def build_summary(slopes, covariance):
    grouped = defaultdict(list)
    for entry in slopes:
        grouped[(entry["eta1"], entry["eta2"])].append(entry)

    summary = []
    for (eta1, eta2), group in sorted(grouped.items()):
        mean, reduced_chi2, spread, count = compatibility(group, covariance)
        summary.append(
            {
                "eta1": eta1,
                "eta2": eta2,
                "n_pt_cells": count,
                "shared_slope": mean,
                "slope_std": spread,
                "reduced_chi2": reduced_chi2,
            }
        )
    return grouped, summary


def make_grid(summary, field):
    neta1 = max(entry["eta1"] for entry in summary) + 1
    neta2 = max(entry["eta2"] for entry in summary) + 1
    grid = np.full((neta1, neta2), np.nan)
    for entry in summary:
        grid[entry["eta1"], entry["eta2"]] = entry[field]
    return grid


def plot_summary(summary, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), constrained_layout=True)
    panels = [
        ("n_pt_cells", "Active pt cells", None),
        ("slope_std", "Slope standard deviation", None),
        ("reduced_chi2", "Compatibility with shared slope: chi2 / ndof", 10.0),
    ]
    for ax, (field, title, vmax) in zip(axes, panels):
        image = ax.imshow(
            make_grid(summary, field).T,
            origin="lower",
            aspect="auto",
            interpolation="nearest",
            vmax=vmax,
        )
        ax.set_xlabel("eta1 bin")
        ax.set_ylabel("eta2 bin")
        ax.set_title(title)
        fig.colorbar(image, ax=ax)
    fig.savefig(outdir / "background_slope_compatibility.png", dpi=160)
    plt.close(fig)


def plot_detail(group, eta1, eta2, outdir):
    npt1 = max(entry["pt1"] for entry in group) + 1
    npt2 = max(entry["pt2"] for entry in group) + 1
    values = np.full((npt1, npt2), np.nan)
    errors = np.full((npt1, npt2), np.nan)
    for entry in group:
        values[entry["pt1"], entry["pt2"]] = entry["value"]
        errors[entry["pt1"], entry["pt2"]] = np.sqrt(entry["variance"])

    fig, ax = plt.subplots(figsize=(7, 5.5), constrained_layout=True)
    image = ax.imshow(values.T, origin="lower", aspect="auto", interpolation="nearest")
    for pt1, pt2 in zip(*np.nonzero(np.isfinite(values))):
        ax.text(
            pt1,
            pt2,
            f"{values[pt1, pt2]:.2g}\n+/- {errors[pt1, pt2]:.2g}",
            ha="center",
            va="center",
            fontsize=8,
        )
    ax.set_xlabel("pt1 bin")
    ax.set_ylabel("pt2 bin")
    ax.set_title(f"Fitted background slopes: eta1={eta1}, eta2={eta2}")
    fig.colorbar(image, ax=ax, label="Exponential slope")
    fig.savefig(outdir / f"background_slopes_eta1_{eta1}_eta2_{eta2}.png", dpi=160)
    plt.close(fig)


def write_summary(summary, outdir):
    fields = ["eta1", "eta2", "n_pt_cells", "shared_slope", "slope_std", "reduced_chi2"]
    with open(
        outdir / "background_slope_compatibility.csv", "w", newline=""
    ) as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)


def main():
    args = parse_args()
    if (args.eta1 is None) != (args.eta2 is None):
        raise ValueError("--eta1 and --eta2 must be passed together")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    slopes, covariance = load_slopes(args.fitresult, args.process, args.useCovariance)
    grouped, summary = build_summary(slopes, covariance)
    ranked = sorted(
        summary,
        key=lambda entry: np.nan_to_num(entry["reduced_chi2"], nan=-np.inf),
        reverse=True,
    )

    write_summary(summary, outdir)
    plot_summary(summary, outdir)
    selected = (args.eta1, args.eta2) if args.eta1 is not None else None
    if selected is None:
        selected = (ranked[0]["eta1"], ranked[0]["eta2"])
    if selected not in grouped:
        raise ValueError(f"No fitted slopes found for eta bins {selected}")
    plot_detail(grouped[selected], *selected, outdir)

    mode = "full covariance" if covariance is not None else "diagonal parameter errors"
    print(f"Found {len(slopes)} fitted slopes; compatibility uses {mode}.")
    print("Largest chi2 / ndof values for sharing slopes across pt bins:")
    for entry in ranked[: args.top]:
        print(
            f"  eta1={entry['eta1']:2d} eta2={entry['eta2']:2d} "
            f"cells={entry['n_pt_cells']:2d} "
            f"slope={entry['shared_slope']: .4g} "
            f"std={entry['slope_std']:.4g} "
            f"chi2/ndof={entry['reduced_chi2']:.4g}"
        )
    print(f"Plots and CSV written to {outdir}")


if __name__ == "__main__":
    main()
