import argparse
import os

import numpy as np

from rabbit import tensorwriter
from wremnants.postprocessing.rabbit_btojpsik_helpers import (
    _reorder_hist_axes,
    _resolve_plot_labels,
    assert_matching_axes,
    collapse_axes,
    load_histogram,
    plot_curvature_response,
    plot_variation_projection,
    rebin_histogram,
    rebin_variation_unc_axis,
)
from wremnants.utilities import common, parsing
from wums import boostHistHelpers as hh


def parse_args():
    analysis_label = common.analysis_label(os.path.basename(__file__))
    parser, initargs = parsing.common_parser(analysis_label)
    parser.description = "Convert the BuToJpsiK histograms into a Rabbit tensor."
    parser.add_argument(
        "-i", "--infile", required=True, help="File that stores the histograms."
    )
    parser.add_argument(
        "--channel", default="btojpsik_stuff", help="Channel name stored in the tensor."
    )
    parser.add_argument(
        "--datasetSignal",
        required=True,
        help="Dataset key for the nominal signal histogram and its variations.",
    )
    parser.add_argument(
        "--background",
        action="append",
        default=[],
        metavar="PROCESS=DATASET",
        help="Background process (repeatable, format: name=dataset_key).",
    )
    parser.add_argument(
        "--signalProcess",
        default="signal",
        help="Process label used for the signal template.",
    )
    parser.add_argument(
        "--signalNormUncertainty",
        type=float,
        default=None,
        help="Optional lnN uncertainty applied to the signal normalization.",
    )
    parser.add_argument(
        "--flatBkgSeed",
        default="dataResidual",
        choices=["dataResidual", "fixedFraction"],
        help="How to seed the artificial flat background normalization in each (eta, pt, charge) bin.",
    )
    parser.add_argument(
        "--flatBkgFraction",
        type=float,
        default=0.40,
        help="Flat-background fraction relative to the signal yield when --flatBkgSeed fixedFraction is used.",
    )

    parser.add_argument(
        "--systematicLabels",
        nargs="+",
        default=("A", "e", "M"),
        help="Order of the labels embedded in the 'unc' axis.",
    )
    parser.add_argument(
        "--variationDownIndex",
        type=int,
        default=0,
        help="Index selecting the down variation.",
    )
    parser.add_argument(
        "--variationUpIndex",
        type=int,
        default=1,
        help="Index selecting the up variation.",
    )

    parser.add_argument(
        "--massAxis", default="bkmm_jpsimc_mass", help="Name of the mass axis."
    )
    parser.add_argument(
        "--ptAxis",
        default="bkmm_kaon_pt",
        help="Name of the kaon-pt axis kept in the tensor.",
    )
    parser.add_argument(
        "--etaAxis",
        default="bkmm_kaon_eta",
        help="Name of the eta axis kept in the tensor.",
    )
    parser.add_argument(
        "--chargeAxis",
        default="bkmm_kaon_charge",
        help="Name of the charge axis kept in the tensor.",
    )
    parser.add_argument(
        "--systematicAxis",
        default="unc",
        help="Axis that enumerates (A,e,M)*neta bins.",
    )
    parser.add_argument(
        "--variationAxis",
        default="downUpVar",
        help="Axis that encodes up/down variations.",
    )
    parser.add_argument(
        "--curvatureAxis",
        default="bkmm_kaon_curvature",
        help="Name of the curvature axis (k) used in the diagnostic plot.",
    )
    parser.add_argument(
        "--systematicType",
        default="log_normal",
        choices=["log_normal", "normal"],
        help="TensorWriter systematic type.",
    )
    parser.add_argument(
        "--plotCurvatureResponse",
        nargs="*",
        metavar="SYST",
        default=None,
        help="Produce k'/k response plots per eta bin. Optionally provide a subset of systematics (e.g. A M) to plot.",
    )
    parser.add_argument(
        "--plotOutput",
        default=None,
        help="Directory where curvature-response plots are stored.",
    )
    parser.add_argument(
        "--plotCurvatureScale",
        type=float,
        default=1.0,
        help="Scale factor applied to (k'/k - 1); useful for magnifying or shrinking variations.",
    )
    parser.add_argument(
        "--plotVariation",
        default=None,
        help="Plot up/down variations projected onto the given axis (e.g. bkmm_jpsimc_mass).",
    )
    parser.add_argument(
        "--plotVariationRatio",
        default=None,
        help="Plot up/down variation ratios projected onto the given axis (e.g. bkmm_jpsimc_curvature).",
    )
    parser.add_argument(
        "--binPlotVariation",
        default=None,
        help="If set, plot the variation axis distribution in each bin of this axis.",
    )
    parser.add_argument(
        "--overlayBinVariations",
        action="store_true",
        help="Overlay each bin of --binPlotVariation on the same plot.",
    )
    parser.add_argument(
        "--plotVariationEta",
        type=int,
        default=None,
        help="Optional eta-bin index for variation plots. If unset, plot all eta bins.",
    )
    parser.add_argument(
        "--plotVariationScale",
        type=float,
        default=1.0,
        help="Scale factor applied to variation curves for yield/ratio plots.",
    )
    parser.add_argument(
        "--debugVariation",
        action="store_true",
        help="Print debug summaries for variation templates.",
    )
    parser.add_argument(
        "--etaBins",
        type=int,
        default=7,
        help="Rebin the eta axis to this many bins before writing. Set to the nominal axis size to disable eta rebinning.",
    )
    parser.add_argument(
        "--massBins",
        type=int,
        default=10,
        help="Rebin the mass axis to this many bins before writing. Set to the nominal axis size to disable mass rebinning.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    # The plotting helpers use snake_case attribute names, while this script
    # follows the repo CLI convention of camelCase options.
    args.charge_axis = args.chargeAxis
    args.curvature_axis = args.curvatureAxis
    args.debug_variation = args.debugVariation
    args.eta_axis = args.etaAxis
    args.overlay_bin_variations = args.overlayBinVariations
    args.plot_curvature_response = args.plotCurvatureResponse
    args.plot_curvature_scale = args.plotCurvatureScale
    args.plot_output = args.plotOutput
    args.plot_variation_scale = args.plotVariationScale
    args.systematic_axis = args.systematicAxis
    args.systematic_labels = tuple(args.systematicLabels)
    args.tensor_output = args.outfolder
    args.variation_axis = args.variationAxis
    args.variation_down_index = args.variationDownIndex
    args.variation_up_index = args.variationUpIndex

    outname = os.path.splitext(os.path.basename(__file__))[0]
    if args.postfix:
        outname += f"_{args.postfix}"

    signal_hist, data_hist, variation_hist = load_histogram(
        args.infile, args.datasetSignal
    )

    print(signal_hist)
    print(data_hist)
    print(variation_hist)

    n_pt_bins = signal_hist.axes[args.ptAxis].size
    rebinning = {}
    if args.etaBins is not None and signal_hist.axes[args.etaAxis].size != args.etaBins:
        rebinning[args.etaAxis] = args.etaBins
    if (
        args.massBins is not None
        and signal_hist.axes[args.massAxis].size != args.massBins
    ):
        rebinning[args.massAxis] = args.massBins
    for axis, bins in rebinning.items():
        print(f"Rebinning {axis} into {bins} bins")

    eta_rebin_factor = None
    if args.etaAxis in rebinning:
        eta_bins = rebinning[args.etaAxis]
        eta_rebin_factor = signal_hist.axes[args.etaAxis].size // eta_bins

    signal_hist = rebin_histogram(signal_hist, rebinning)
    data_hist = rebin_histogram(data_hist, rebinning)

    variation_rebinning = dict(rebinning)
    if eta_rebin_factor is not None and args.systematicAxis in variation_hist.axes.name:
        variation_hist = rebin_variation_unc_axis(
            variation_hist,
            args.etaAxis,
            args.systematicAxis,
            args.systematicLabels,
            eta_rebin_factor,
        )
        variation_rebinning.pop(args.etaAxis, None)
    variation_hist = rebin_histogram(variation_hist, variation_rebinning)

    # lose stats quick at "high" pT...

    # import pdb
    # pdb.set_trace()

    background_hists = {}
    for spec in args.background:
        if "=" not in spec:
            raise argparse.ArgumentTypeError(
                f"Background specification '{spec}' does not contain '='."
            )
        process, dataset = spec.split("=", maxsplit=1)
        background_nominal, _, _ = load_histogram(args.infile, dataset)
        background_hists[process] = background_nominal

    # drop_nominal_axes = [args.ptAxis, args.systematicAxis, args.variationAxis]
    # drop_nominal_axes = [args.systematicAxis, args.variationAxis]
    # signal_hist = collapse_axes(signal_hist, drop_nominal_axes)
    # for key, h in background_hists.items():
    #    background_hists[key] = collapse_axes(h, drop_nominal_axes)

    assert_matching_axes(background_hists, signal_hist, label="Background")

    n_eta_bins = signal_hist.axes[args.etaAxis].size
    n_charge_bins = signal_hist.axes[args.chargeAxis].size

    bkg_hist = signal_hist.copy()
    bkg_hist.values()[...] = 0.0
    if bkg_hist.variances() is not None:
        bkg_hist.variances()[...] = 0.0

    mass_axis_idx = signal_hist.axes.name.index(args.massAxis)
    mass_bins = signal_hist.axes[args.massAxis].size
    floor_yield = 1e-6

    for icharge in range(n_charge_bins):
        for ipt in range(n_pt_bins):
            for ieta in range(n_eta_bins):
                sel = {
                    args.chargeAxis: icharge,
                    args.ptAxis: ipt,
                    args.etaAxis: ieta,
                }
                signal_yield = signal_hist[sel].values().sum()
                data_yield = data_hist[sel].values().sum()

                if args.flatBkgSeed == "dataResidual":
                    seed_yield = max(data_yield - signal_yield, 0.0)
                else:
                    seed_yield = args.flatBkgFraction * signal_yield

                if seed_yield <= 0.0:
                    continue

                target = [slice(None)] * bkg_hist.values().ndim
                target[signal_hist.axes.name.index(args.chargeAxis)] = icharge
                target[signal_hist.axes.name.index(args.ptAxis)] = ipt
                target[signal_hist.axes.name.index(args.etaAxis)] = ieta
                target[mass_axis_idx] = slice(None)
                target = tuple(target)

                bkg_hist.values()[target] = max(seed_yield / mass_bins, floor_yield)

    # tensor writer now
    writer = tensorwriter.TensorWriter(systematic_type=args.systematicType)
    writer.add_channel(signal_hist.axes, name=args.channel)
    writer.add_data(data_hist, channel=args.channel)
    # This calibration tensor already floats the signal yield independently in every
    # fitted (eta, pt, charge) bin below. Marking the same process as a Rabbit
    # "signal" would add a redundant global signal-strength POI on top of those
    # bin-wise yields, which makes the explicit-fit Hessian singular.
    writer.add_process(signal_hist, args.signalProcess, args.channel, signal=False)
    for proc, h in background_hists.items():
        writer.add_process(h, proc, args.channel)
    # artificial background
    writer.add_process(bkg_hist, "flatBkg", args.channel)

    # if args.signalNormUncertainty is not None:
    #    writer.add_norm_systematic(
    #        name="signal_norm",
    #        process=args.signalProcess,
    #        channel=args.channel,
    #        uncertainty=args.signalNormUncertainty,
    #    )

    n_eta_bins = signal_hist.axes[args.etaAxis].size
    n_charge_bins = signal_hist.axes[args.chargeAxis].size

    # free float yields in bins of pt eta charge
    new_sig_basis = hh.expand_hist_by_duplicate_axes(
        signal_hist,
        [args.chargeAxis, args.ptAxis, args.etaAxis],
        ["chargeVar", "ptVar", "etaVar"],
        put_trailing=True,
        flow=False,
    )
    new_bkg_basis = hh.expand_hist_by_duplicate_axes(
        bkg_hist,
        [args.chargeAxis, args.ptAxis, args.etaAxis],
        ["chargeVar", "ptVar", "etaVar"],
        put_trailing=True,
        flow=False,
    )

    procs = {args.signalProcess: signal_hist, "flatBkg": bkg_hist}
    basis_by_proc = {
        args.signalProcess: new_sig_basis,
        "flatBkg": new_bkg_basis,
    }
    unc = 0.1
    for proc, nominal_hist in procs.items():
        basis_hist = basis_by_proc[proc]
        for icharge in range(n_charge_bins):
            for ipt in range(n_pt_bins):
                for ieta in range(n_eta_bins):

                    mask = {"etaVar": ieta, "ptVar": ipt, "chargeVar": icharge}
                    bin = basis_hist[mask]
                    if not np.any(bin.values()):
                        continue

                    syst_name = f"norm_{proc}_eta{ieta}_pt{ipt}_charge{icharge}"

                    up_hist = nominal_hist + bin * unc
                    # down_hist--symmetric

                    writer.add_systematic(
                        up_hist,
                        name=syst_name,
                        process=proc,
                        channel=args.channel,
                        constrained=False,
                        noi=True,
                    )

    # signal yield systematics matching binning for A,e,M
    # for ieta in range(n_eta_bins):
    #    mask = {"etaVar": ieta}
    #    bin = new_sig_basis[mask]
    #    syst_name = f"norm_{args.signalProcess}_eta{ieta}"
    #    bin = new_sig_basis[{"etaVar": ieta}]
    #    up_hist = signal_hist + bin * unc
    #
    #    writer.add_systematic(
    #        up_hist,
    #        name=syst_name,
    #        process=args.signalProcess,
    #        channel=args.channel,
    #        constrained=False,
    #        noi=True
    #    )

    # A e M bitchhhh
    labels = tuple(args.systematicLabels)

    plot_labels = _resolve_plot_labels(args)
    if plot_labels:
        plot_curvature_response(signal_hist, variation_hist, args, plot_labels)
    if args.plotVariation and args.plotVariationRatio:
        raise ValueError(
            "--plotVariation and --plotVariationRatio are mutually exclusive."
        )
    if args.plotVariation:
        plot_variation_projection(
            signal_hist,
            variation_hist,
            args,
            args.plotVariation,
            eta_only=args.plotVariationEta,
            mode="yield",
            bin_axis=args.binPlotVariation,
        )
    if args.plotVariationRatio:
        plot_variation_projection(
            signal_hist,
            variation_hist,
            args,
            args.plotVariationRatio,
            eta_only=args.plotVariationEta,
            mode="ratio",
            bin_axis=args.binPlotVariation,
        )

    for eta_idx in range(n_eta_bins):
        for label in labels:
            offset = labels.index(label)
            systematic_bin = len(labels) * eta_idx + offset

            up_selection = {
                args.systematicAxis: systematic_bin,
                args.variationAxis: args.variationUpIndex,
            }
            down_selection = {
                args.systematicAxis: systematic_bin,
                args.variationAxis: args.variationDownIndex,
            }

            up_variation = collapse_axes(
                variation_hist[up_selection],
                [args.systematicAxis, args.variationAxis],
            )
            down_variation = collapse_axes(
                variation_hist[down_selection],
                [args.systematicAxis, args.variationAxis],
            )

            target_axes = signal_hist.axes.name
            up_variation = _reorder_hist_axes(up_variation, target_axes)
            down_variation = _reorder_hist_axes(down_variation, target_axes)

            if up_variation.axes.name != signal_hist.axes.name:
                raise RuntimeError(
                    f"Up variation axes {up_variation.axes.name} do not match nominal axes {signal_hist.axes.name}."
                )
            if down_variation.axes.name != signal_hist.axes.name:
                raise RuntimeError(
                    f"Down variation axes {down_variation.axes.name} do not match nominal axes {signal_hist.axes.name}."
                )

            up_hist = signal_hist.copy()
            down_hist = signal_hist.copy()

            # Each A/e/M nuisance is defined for one eta bin. The selected
            # variation slice only carries that eta bin, so keep the nominal
            # template elsewhere and replace only the targeted eta slice.
            eta_target = [slice(None)] * up_hist.values().ndim
            eta_axis_idx = up_hist.axes.name.index(args.etaAxis)
            eta_target[eta_axis_idx] = eta_idx
            eta_target = tuple(eta_target)

            up_hist.values()[eta_target] = up_variation.values()[eta_target]
            down_hist.values()[eta_target] = down_variation.values()[eta_target]

            up_vars = up_variation.variances()
            down_vars = down_variation.variances()
            if up_vars is not None:
                up_hist.variances()[eta_target] = up_vars[eta_target]
            if down_vars is not None:
                down_hist.variances()[eta_target] = down_vars[eta_target]

            writer.add_systematic(
                [up_hist, down_hist],
                name=f"{label}_eta{eta_idx}",
                process=args.signalProcess,
                channel=args.channel,
                constrained=False,
                noi=True,
            )
        # up_hist.values()[..., charge_idx] = up_variation.values()[..., charge_idx]
        # down_hist.values()[..., charge_idx] = down_variation.values()[..., charge_idx]
    #
    # up_vars = up_variation.variances()
    # down_vars = down_variation.variances()
    # up_hist_vars = up_hist.variances()
    # down_hist_vars = down_hist.variances()
    #
    # if up_vars is not None and up_hist_vars is not None:
    #    up_hist_vars[..., charge_idx] = up_vars[..., charge_idx]
    # if down_vars is not None and down_hist_vars is not None:
    #    down_hist_vars[..., charge_idx] = down_vars[..., charge_idx]
    #
    # writer.add_systematic(
    #    [up_hist, down_hist],
    #    name=f"{label}_eta{eta_idx}",
    #    process=args.signalProcess,
    #    channel=args.channel,
    #    constrained=False,
    # )

    writer.write(args.outfolder, outname)


if __name__ == "__main__":
    main()
