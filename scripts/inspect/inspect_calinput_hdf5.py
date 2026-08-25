#!/usr/bin/env python3

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from wremnants.utilities.io_tools import base_io
from wums import ioutils

SCALE_UNC_HIST = "nominal_muonScaleSyst_responseWeights"
SCALE_UNC_LABELS = ("A", "e", "M")
SCALE_DOWN_UP_LABELS = ("Down", "Up")
ETA_REGION_BINS = {
    "barrel": 12,
    "endcap": 23,
}


def materialize(obj):
    return obj.get() if isinstance(obj, ioutils.H5PickleProxy) else obj


def hist_sum(h):
    values = h.values(flow=False)
    return float(np.sum(values))


def print_hist_summary(name, h):
    print(f"    hist: {name}")
    print(f"      axes: {[axis.name for axis in h.axes]}")
    print(f"      shape: {h.values(flow=False).shape}")
    print(f"      sum: {hist_sum(h):.6g}")
    for axis in h.axes:
        edges = axis.edges
        print(
            f"      axis {axis.name}: bins={axis.size}, "
            f"range=[{edges[0]:.6g}, {edges[-1]:.6g}]"
        )


def load_results(path):
    h5file = h5py.File(path, "r")
    results = base_io.load_results_h5py(h5file)

    # Force lazy histogram proxies to load while the HDF5 file is still open.
    for result in results.values():
        if not isinstance(result, dict) or "output" not in result:
            continue
        for name, obj in list(result["output"].items()):
            result["output"][name] = materialize(obj)

    h5file.close()
    return results


def get_histograms(results):
    histograms = {}
    for dataset_name, result in results.items():
        if dataset_name == "meta_info" or not isinstance(result, dict):
            continue
        for hist_name, hist_obj in result.get("output", {}).items():
            histograms[(dataset_name, hist_name)] = hist_obj
    return histograms


def step_plot(ax, h1d, label):
    axis = h1d.axes[0]
    edges = axis.edges
    values = h1d.values(flow=False)
    step_plot_values(ax, edges, values, axis.name, label)


def step_plot_values(ax, edges, values, xlabel, label):
    y = np.r_[values, values[-1] if len(values) else 0.0]
    ax.step(edges, y, where="post", label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Events")


def project_values(h, keep_names, selections=None):
    names = list(h.axes.name)
    selections = selections or {}
    slicer = []
    remaining_names = []

    for name in names:
        if name in selections:
            slicer.append(selections[name])
        else:
            slicer.append(slice(None))
            remaining_names.append(name)

    values = np.asarray(h.values(flow=False)[tuple(slicer)])

    for axis_index in reversed(range(values.ndim)):
        if remaining_names[axis_index] not in keep_names:
            values = values.sum(axis=axis_index)
            del remaining_names[axis_index]

    if remaining_names != list(keep_names):
        order = [remaining_names.index(name) for name in keep_names]
        values = np.moveaxis(values, order, range(len(order)))

    return values


def plot_projection(histograms, axis_name, outpath):
    fig, ax = plt.subplots(figsize=(7, 5))
    drew = False

    for (_, hist_name), h in sorted(histograms.items()):
        if hist_name == SCALE_UNC_HIST:
            continue
        if axis_name not in h.axes.name:
            continue
        values = project_values(h, [axis_name])
        step_plot_values(ax, h.axes[axis_name].edges, values, axis_name, hist_name)
        drew = True

    if not drew:
        plt.close(fig)
        return

    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def plot_2d(histograms, xaxis, yaxis, outdir):
    for (_, hist_name), h in sorted(histograms.items()):
        if hist_name == SCALE_UNC_HIST:
            continue
        if xaxis not in h.axes.name or yaxis not in h.axes.name:
            continue
        xedges = h.axes[xaxis].edges
        yedges = h.axes[yaxis].edges
        values = project_values(h, [xaxis, yaxis]).T

        fig, ax = plt.subplots(figsize=(6, 5))
        mesh = ax.pcolormesh(xedges, yedges, values)
        ax.set_xlabel(xaxis)
        ax.set_ylabel(yaxis)
        ax.set_title(hist_name)
        fig.colorbar(mesh, ax=ax, label="Events")
        fig.tight_layout()
        fig.savefig(outdir / f"{hist_name}_{xaxis}_{yaxis}.png")
        plt.close(fig)


def nominal_hist_for_systematic(histograms, dataset_name):
    for (candidate_dataset, hist_name), h in histograms.items():
        if candidate_dataset != dataset_name:
            continue
        if hist_name == SCALE_UNC_HIST:
            continue
        if {"eta1", "eta2", "pt1", "pt2", "mass"}.issubset(set(h.axes.name)):
            return h
    return None


def selected_eta_parameter_unc_index(hsyst, param_index, eta_bin):
    n_unc = hsyst.axes["unc"].size
    unc_index = eta_bin * len(SCALE_UNC_LABELS) + param_index
    if unc_index >= n_unc:
        return None
    return unc_index


def plot_muon_scale_systematics(histograms, outdir):
    for (dataset_name, hist_name), hsyst in sorted(histograms.items()):
        if hist_name != SCALE_UNC_HIST:
            continue
        if "unc" not in hsyst.axes.name or "downUpVar" not in hsyst.axes.name:
            print(f"  Skipping {hist_name}: expected tensor axes 'unc' and 'downUpVar'")
            continue

        hnom = nominal_hist_for_systematic(histograms, dataset_name)
        if hnom is None:
            print(f"  Skipping {hist_name}: no matching nominal histogram found")
            continue

        nom_values = project_values(hnom, ["mass"])
        mass_edges = hnom.axes["mass"].edges
        xcenters = 0.5 * (mass_edges[:-1] + mass_edges[1:])
        variation_cache = {}
        n_eta_unc_bins = hsyst.axes["unc"].size // len(SCALE_UNC_LABELS)
        if hsyst.axes["unc"].size % len(SCALE_UNC_LABELS) != 0:
            print(
                f"  Warning: {hist_name} has {hsyst.axes['unc'].size} unc bins, "
                "which is not divisible by 3. Eta-bin A/e/M selection assumes calinput "
                "prefit output made with --fitMuonScaleAndResolution."
            )

        for region, requested_eta_bin in ETA_REGION_BINS.items():
            eta_bin = min(requested_eta_bin, max(n_eta_unc_bins - 1, 0))
            if eta_bin != requested_eta_bin:
                print(
                    f"  Warning: requested {region} eta-bin {requested_eta_bin}, "
                    f"but only {n_eta_unc_bins} eta bins are available. Using {eta_bin}."
                )

            for iunc, unc_label in enumerate(SCALE_UNC_LABELS):
                unc_index = selected_eta_parameter_unc_index(hsyst, iunc, eta_bin)
                if unc_index is None:
                    continue

                for idownup in range(len(SCALE_DOWN_UP_LABELS)):
                    variation_cache[(unc_index, idownup)] = project_values(
                        hsyst,
                        ["mass"],
                        {"unc": unc_index, "downUpVar": idownup},
                    )

                fig, (ax, rax) = plt.subplots(
                    2,
                    1,
                    figsize=(7, 6),
                    sharex=True,
                    gridspec_kw={"height_ratios": [3, 1]},
                )
                step_plot_values(ax, mass_edges, nom_values, "mass", "nominal")

                for idownup, direction in enumerate(SCALE_DOWN_UP_LABELS):
                    values = variation_cache[(unc_index, idownup)]
                    y = np.r_[values, values[-1] if len(values) else 0.0]
                    ax.step(
                        mass_edges,
                        y,
                        where="post",
                        label=f"{unc_label} {direction}",
                    )

                    ratio = np.divide(
                        values,
                        nom_values,
                        out=np.ones_like(values, dtype=float),
                        where=nom_values != 0,
                    )
                    rax.plot(
                        xcenters, ratio, marker="o", linestyle="-", label=direction
                    )

                ax.set_ylabel("Events")
                ax.set_title(
                    f"{dataset_name}: muon scale {unc_label}, {region} eta bin {eta_bin}"
                )
                ax.legend()
                ax.grid(True, alpha=0.25)
                rax.axhline(1.0, color="black", linewidth=1)
                rax.set_xlabel("mass")
                rax.set_ylabel("var/nom")
                rax.grid(True, alpha=0.25)
                fig.tight_layout()
                fig.savefig(
                    outdir
                    / f"{dataset_name}_muonScale_{unc_label}_{region}_etaBin{eta_bin}_mass.png"
                )
                plt.close(fig)

                for idownup, direction in enumerate(SCALE_DOWN_UP_LABELS):
                    values = variation_cache[(unc_index, idownup)]
                    ratio = np.divide(
                        values,
                        nom_values,
                        out=np.ones_like(values, dtype=float),
                        where=nom_values != 0,
                    )

                    fig, (ax, rax) = plt.subplots(
                        2,
                        1,
                        figsize=(7, 6),
                        sharex=True,
                        gridspec_kw={"height_ratios": [3, 1]},
                    )
                    step_plot_values(ax, mass_edges, nom_values, "mass", "nominal")
                    y = np.r_[values, values[-1] if len(values) else 0.0]
                    ax.step(
                        mass_edges,
                        y,
                        where="post",
                        label=f"{unc_label} {direction}",
                    )
                    ax.set_ylabel("Events")
                    ax.set_title(
                        f"{dataset_name}: muon scale {unc_label} {direction}, "
                        f"{region} eta bin {eta_bin}"
                    )
                    ax.legend()
                    ax.grid(True, alpha=0.25)

                    rax.plot(xcenters, ratio, marker="o", linestyle="-")
                    rax.axhline(1.0, color="black", linewidth=1)
                    rax.set_xlabel("mass")
                    rax.set_ylabel("var/nom")
                    rax.grid(True, alpha=0.25)
                    fig.tight_layout()
                    fig.savefig(
                        outdir
                        / f"{dataset_name}_muonScale_{unc_label}_{direction}_{region}_etaBin{eta_bin}_mass.png"
                    )
                    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Inspect calInput-style dimuon resonance HDF5 output and make simple plots."
    )
    parser.add_argument(
        "input", help="Input HDF5 file from dimuon_resonances_calinput.py"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default=None,
        help="Directory for plots. Defaults to '<input stem>_plots'.",
    )
    args = parser.parse_args()

    infile = Path(args.input)
    outdir = (
        Path(args.outdir)
        if args.outdir
        else infile.with_suffix("").parent / (infile.with_suffix("").name + "_plots")
    )
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading {infile}")
    results = load_results(infile)

    print("\nDatasets")
    for dataset_name, result in results.items():
        if dataset_name == "meta_info" or not isinstance(result, dict):
            continue
        dataset = result.get("dataset", {})
        print(f"  {dataset_name}")
        print(f"    is_data: {dataset.get('is_data')}")
        print(f"    files: {dataset.get('filepaths')}")
        print(f"    event_count: {result.get('event_count')}")
        print(f"    weight_sum: {result.get('weight_sum')}")

    histograms = get_histograms(results)
    print("\nHistograms")
    for (_, hist_name), h in sorted(histograms.items()):
        print_hist_summary(hist_name, h)

    for axis_name in ["mass", "pt1", "pt2", "eta1", "eta2"]:
        plot_projection(histograms, axis_name, outdir / f"{axis_name}.png")

    plot_2d(histograms, "eta1", "eta2", outdir)
    plot_2d(histograms, "pt1", "pt2", outdir)
    plot_muon_scale_systematics(histograms, outdir)

    print(f"\nWrote plots to {outdir}")


if __name__ == "__main__":
    main()
