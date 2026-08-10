import argparse
import os

import numpy as np

from rabbit import tensorwriter
from wremnants.postprocessing.rabbit_btojpsik_helpers import (
    _reorder_hist_axes,
    _resolve_plot_labels,
    assert_matching_axes,
    collapse_axes,
    evaluate_exp_plus_constant,
    fit_exp_plus_constant_seed,
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
        help="How to seed the artificial background normalization in each (eta, pt, charge) bin.",
    )
    parser.add_argument(
        "--flatBkgFraction",
        type=float,
        default=0.40,
        help="Flat-background fraction relative to the signal yield when --flatBkgSeed fixedFraction is used.",
    )
    parser.add_argument(
        "--bkgModel",
        default="flat",
        choices=["exp", "flat", "exp_poi"],
        help="Background model. 'flat': legacy constant background with per-bin NOIs. "
        "'exp': A*exp(-B*x)+C with per-bin NOIs (deprecated). "
        "'exp_poi': single all-ones bkgExp process; shape exp(lnAmpl+slope*x) "
        "handled entirely by AxisExpModel POIs (no flatBkg process).",
    )
    parser.add_argument(
        "--signalNormPOI",
        action="store_true",
        help="Signal per-bin yields are handled by AxisNormModel (or another "
        "POI model) instead of template-morph NOIs. When set, per-bin signal "
        "norm systematics are NOT added to the tensor and the signal process "
        "is marked signal=True; A/e/M systematics are still added as NOIs.",
    )
    parser.add_argument(
        "--bkgPOI",
        action="store_true",
        help="Write a single all-ones 'background' process whose shape and "
        "amplitude are fully owned by a POI model (e.g. AxisExpModel or "
        "AxisBernsteinModel). No background NOIs are added to the tensor.",
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
        "--nomc",
        action="store_true",
        help="Switch to nomc mode: use bkmm_nomc_mass as fit observable and add "
        "A/e/M nuisances for each J/psi muon.",
    )
    parser.add_argument(
        "--massAxis",
        default=None,
        help="Name of the mass axis. Defaults to bkmm_nomc_mass when --nomc is set, "
        "bkmm_jpsimc_mass otherwise.",
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
        "--normalizeAeM",
        action="store_true",
        help="Normalize A/e/M variation templates so the mass-integrated yield per "
        "(pt, eta, charge) cell equals the nominal. Makes A/e/M pure mass-shape "
        "variations, removing the normalization redundancy with AxisNormModel POIs "
        "and preventing non-positive-definite Hessian failures when using "
        "--signalNormPOI.",
    )
    parser.add_argument(
        "--etaBins",
        type=int,
        default=None,
        help="Rebin the eta axis to this many bins before writing. Set to the nominal axis size to disable eta rebinning.",
    )
    parser.add_argument(
        "--massBins",
        type=int,
        default=10,
        help="Rebin the mass axis to this many bins before writing. Set to the nominal axis size to disable mass rebinning.",
    )
    return parser.parse_args()


def _add_aem_systematics(
    writer, variation_hist, signal_hist, args, labels, n_eta_bins, name_suffix=""
):
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
                    f"Up variation axes {up_variation.axes.name} do not match "
                    f"nominal axes {signal_hist.axes.name}."
                )
            if down_variation.axes.name != signal_hist.axes.name:
                raise RuntimeError(
                    f"Down variation axes {down_variation.axes.name} do not match "
                    f"nominal axes {signal_hist.axes.name}."
                )

            up_hist = signal_hist.copy()
            down_hist = signal_hist.copy()

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

            if args.normalizeAeM:
                mass_ax = up_hist.axes.name.index(args.massAxis)
                # After eta_target fixes eta with an integer index, the mass axis
                # position in the reduced array shifts left by 1 if it was after eta.
                mass_ax_r = mass_ax if mass_ax < eta_axis_idx else mass_ax - 1
                nom_sum = signal_hist.values()[eta_target].sum(
                    axis=mass_ax_r, keepdims=True
                )
                up_slice = up_hist.values()[eta_target].copy()
                up_sum = up_slice.sum(axis=mass_ax_r, keepdims=True)
                up_scale = np.where(up_sum > 0, nom_sum / up_sum, 1.0)
                up_hist.values()[eta_target] = up_slice * up_scale

                down_slice = down_hist.values()[eta_target].copy()
                down_sum = down_slice.sum(axis=mass_ax_r, keepdims=True)
                down_scale = np.where(down_sum > 0, nom_sum / down_sum, 1.0)
                down_hist.values()[eta_target] = down_slice * down_scale

            writer.add_systematic(
                [up_hist, down_hist],
                name=f"{label}{name_suffix}_eta{eta_idx}",
                process=args.signalProcess,
                channel=args.channel,
                constrained=False,
                noi=True,
            )


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
    if args.massAxis is None:
        args.massAxis = "bkmm_nomc_mass" if args.nomc else "bkmm_jpsimc_mass"

    outname = os.path.splitext(os.path.basename(__file__))[0]
    if args.postfix:
        outname += f"_{args.postfix}"

    signal_hist, data_hist, variation_hist = load_histogram(
        args.infile, args.datasetSignal
    )

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

    bkg_exp_hist = None
    bkg_const_hist = None
    if args.bkgModel == "exp":
        bkg_exp_hist = signal_hist.copy()
        bkg_exp_hist.values()[...] = 0.0
        if bkg_exp_hist.variances() is not None:
            bkg_exp_hist.variances()[...] = 0.0

        bkg_const_hist = signal_hist.copy()
        bkg_const_hist.values()[...] = 0.0
        if bkg_const_hist.variances() is not None:
            bkg_const_hist.variances()[...] = 0.0

    # exp_poi mode: single all-ones template; shape owned by AxisExpModel POIs
    bkg_exp_poi_hist = None
    if args.bkgModel == "exp_poi":
        bkg_exp_poi_hist = signal_hist.copy()
        bkg_exp_poi_hist.values()[...] = 1.0
        if bkg_exp_poi_hist.variances() is not None:
            bkg_exp_poi_hist.variances()[...] = 0.0

    # --bkgPOI: generic all-ones template; shape owned by AxisExpModel/AxisBernsteinModel
    bkg_poi_hist = None
    if args.bkgPOI:
        bkg_poi_hist = signal_hist.copy()
        bkg_poi_hist.values()[...] = 1.0
        if bkg_poi_hist.variances() is not None:
            bkg_poi_hist.variances()[...] = 0.0

    mass_axis_idx = signal_hist.axes.name.index(args.massAxis)
    mass_bins = signal_hist.axes[args.massAxis].size
    mass_coordinate = np.asarray(signal_hist.axes[args.massAxis].centers, dtype=float)
    mass_coordinate = mass_coordinate - mass_coordinate[0]
    mass_span = max(float(mass_coordinate[-1]), 1e-6)
    floor_yield = 1e-6
    # These background coefficients are free fit parameters, so use small local
    # finite-difference steps to define their tensor directions rather than a
    # full-component jump as the nominal prefit variation.
    bkg_amp_step_fraction = 0.10

    bkg_seed_A = np.zeros((n_charge_bins, n_pt_bins, n_eta_bins), dtype=float)
    bkg_seed_C = np.zeros((n_charge_bins, n_pt_bins, n_eta_bins), dtype=float)
    bkg_seed_total = np.zeros((n_charge_bins, n_pt_bins, n_eta_bins), dtype=float)

    background_process = "flatBkg" if args.bkgModel == "flat" else "expBkg"

    for icharge in range(n_charge_bins):
        for ipt in range(n_pt_bins):
            for ieta in range(n_eta_bins):
                sel = {
                    args.chargeAxis: icharge,
                    args.ptAxis: ipt,
                    args.etaAxis: ieta,
                }
                signal_slice = signal_hist[sel]
                data_slice = data_hist[sel]
                signal_yield = signal_slice.values().sum()
                data_yield = data_slice.values().sum()
                residual_spectrum = np.clip(
                    data_slice.values() - signal_slice.values(), 0.0, None
                )

                if args.flatBkgSeed == "dataResidual":
                    seed_yield = float(np.sum(residual_spectrum))
                else:
                    seed_yield = args.flatBkgFraction * signal_yield
                    residual_spectrum = None

                if seed_yield <= 0.0:
                    continue

                target = [slice(None)] * bkg_hist.values().ndim
                target[signal_hist.axes.name.index(args.chargeAxis)] = icharge
                target[signal_hist.axes.name.index(args.ptAxis)] = ipt
                target[signal_hist.axes.name.index(args.etaAxis)] = ieta
                target[mass_axis_idx] = slice(None)
                target = tuple(target)

                if args.bkgModel == "flat":
                    model_values = np.full(
                        mass_bins, seed_yield / mass_bins, dtype=float
                    )
                    seed_A = 0.0
                    seed_B = 0.0
                    seed_C = float(seed_yield / mass_bins)
                elif args.bkgModel == "exp_poi" or args.bkgPOI:
                    # All-ones templates; shape is owned by POI models.
                    continue
                else:
                    if (
                        residual_spectrum is not None
                        and np.sum(residual_spectrum) > 0.0
                    ):
                        fit = fit_exp_plus_constant_seed(
                            residual_spectrum, mass_coordinate, floor=floor_yield
                        )
                        seed_A = float(fit["A"])
                        seed_B = max(float(fit["B"]), 0.5 / mass_span)
                        seed_C = float(fit["C"])
                    else:
                        seed_B = 0.5 / mass_span
                        seed_A = 0.7 * seed_yield
                        seed_C = 0.3 * seed_yield / mass_bins
                    component_floor = max(0.1 * seed_yield / mass_bins, floor_yield)
                    exp_component = evaluate_exp_plus_constant(
                        mass_coordinate,
                        max(seed_A, component_floor),
                        seed_B,
                        0.0,
                        floor=floor_yield,
                    )
                    const_component = np.full(
                        mass_bins, max(seed_C, component_floor), dtype=float
                    )
                    component_scale = seed_yield / (
                        np.sum(exp_component) + np.sum(const_component)
                    )
                    exp_component = np.clip(
                        exp_component * component_scale, floor_yield, None
                    )
                    const_component = np.clip(
                        const_component * component_scale, floor_yield, None
                    )
                    model_values = exp_component + const_component
                    bkg_exp_hist.values()[target] = exp_component
                    bkg_const_hist.values()[target] = const_component
                    seed_A = float(max(seed_A, component_floor) * component_scale)
                    seed_C = float(max(seed_C, component_floor) * component_scale)

                model_values = np.clip(model_values, floor_yield, None)
                bkg_hist.values()[target] = model_values
                bkg_seed_A[icharge, ipt, ieta] = seed_A
                bkg_seed_C[icharge, ipt, ieta] = seed_C
                bkg_seed_total[icharge, ipt, ieta] = float(np.sum(model_values))

    # tensor writer now
    writer = tensorwriter.TensorWriter(systematic_type=args.systematicType)
    writer.add_channel(signal_hist.axes, name=args.channel)
    writer.add_data(data_hist, channel=args.channel)
    # When --signalNormPOI is set, per-bin yields are handled by AxisNormModel
    # and no template-morph NOIs are added below, so signal=True is safe.
    # Without --signalNormPOI the legacy NOI path is used and signal=False
    # prevents a redundant global signal-strength POI.
    writer.add_process(
        signal_hist, args.signalProcess, args.channel, signal=args.signalNormPOI
    )
    for proc, h in background_hists.items():
        writer.add_process(h, proc, args.channel)
    # artificial background
    if args.bkgPOI:
        writer.add_process(bkg_poi_hist, "background", args.channel)
    elif args.bkgModel == "flat":
        writer.add_process(bkg_hist, background_process, args.channel)
    elif args.bkgModel == "exp_poi":
        writer.add_process(bkg_exp_poi_hist, "bkgExp", args.channel)
    else:
        # exp model: separate exp/const processes with per-bin systematics
        writer.add_process(bkg_exp_hist, "bkgExp", args.channel)
        writer.add_process(bkg_const_hist, "bkgConst", args.channel)

    # if args.signalNormUncertainty is not None:
    #    writer.add_norm_systematic(
    #        name="signal_norm",
    #        process=args.signalProcess,
    #        channel=args.channel,
    #        uncertainty=args.signalNormUncertainty,
    #    )

    n_eta_bins = signal_hist.axes[args.etaAxis].size
    n_charge_bins = signal_hist.axes[args.chargeAxis].size

    # Free-float yields in bins of pt eta charge.
    # When --signalNormPOI is set, per-bin norms for both signal and flat
    # background are handled by AxisNormModel and NOT added as template-morph
    # systematics here.
    procs = {}
    basis_by_proc = {}
    if not args.signalNormPOI:
        new_sig_basis = hh.expand_hist_by_duplicate_axes(
            signal_hist,
            [args.chargeAxis, args.ptAxis, args.etaAxis],
            ["chargeVar", "ptVar", "etaVar"],
            put_trailing=True,
            flow=False,
        )
        procs[args.signalProcess] = signal_hist
        basis_by_proc[args.signalProcess] = new_sig_basis
    if args.bkgModel == "flat" and not args.signalNormPOI and not args.bkgPOI:
        new_bkg_basis = hh.expand_hist_by_duplicate_axes(
            bkg_hist,
            [args.chargeAxis, args.ptAxis, args.etaAxis],
            ["chargeVar", "ptVar", "etaVar"],
            put_trailing=True,
            flow=False,
        )
        procs[background_process] = bkg_hist
        basis_by_proc[background_process] = new_bkg_basis
    elif args.bkgModel == "exp":
        new_bkg_exp_basis = hh.expand_hist_by_duplicate_axes(
            bkg_exp_hist,
            [args.chargeAxis, args.ptAxis, args.etaAxis],
            ["chargeVar", "ptVar", "etaVar"],
            put_trailing=True,
            flow=False,
        )
        new_bkg_const_basis = hh.expand_hist_by_duplicate_axes(
            bkg_const_hist,
            [args.chargeAxis, args.ptAxis, args.etaAxis],
            ["chargeVar", "ptVar", "etaVar"],
            put_trailing=True,
            flow=False,
        )
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

    if args.bkgModel == "exp":
        charge_axis_idx = bkg_hist.axes.name.index(args.chargeAxis)
        pt_axis_idx = bkg_hist.axes.name.index(args.ptAxis)
        eta_axis_idx = bkg_hist.axes.name.index(args.etaAxis)

        def make_component_variation(hist_template, values, target):
            variation_hist = hist_template.copy()
            variation_hist.values()[target] = np.clip(values, floor_yield, None)
            if variation_hist.variances() is not None:
                variation_hist.variances()[target] = 0.0
            return variation_hist

        for icharge in range(n_charge_bins):
            for ipt in range(n_pt_bins):
                for ieta in range(n_eta_bins):
                    seed_yield = bkg_seed_total[icharge, ipt, ieta]
                    if seed_yield <= 0.0:
                        continue

                    target = [slice(None)] * bkg_hist.values().ndim
                    target[charge_axis_idx] = icharge
                    target[pt_axis_idx] = ipt
                    target[eta_axis_idx] = ieta
                    target[mass_axis_idx] = slice(None)
                    target = tuple(target)

                    seed_A = bkg_seed_A[icharge, ipt, ieta]
                    seed_C = bkg_seed_C[icharge, ipt, ieta]
                    exp_nominal = np.asarray(bkg_exp_hist.values()[target], dtype=float)
                    const_nominal = np.asarray(
                        bkg_const_hist.values()[target], dtype=float
                    )

                    writer.add_systematic(
                        [
                            make_component_variation(
                                bkg_exp_hist,
                                np.clip(
                                    exp_nominal * (1.0 + bkg_amp_step_fraction),
                                    floor_yield,
                                    None,
                                ),
                                target,
                            ),
                            make_component_variation(
                                bkg_exp_hist,
                                np.clip(
                                    exp_nominal * (1.0 - bkg_amp_step_fraction),
                                    floor_yield,
                                    None,
                                ),
                                target,
                            ),
                        ],
                        name=f"bkgA_eta{ieta}_pt{ipt}_charge{icharge}",
                        process="bkgExp",
                        channel=args.channel,
                        constrained=False,
                        noi=True,
                    )
                    writer.add_systematic(
                        [
                            make_component_variation(
                                bkg_const_hist,
                                np.clip(
                                    const_nominal * (1.0 + bkg_amp_step_fraction),
                                    floor_yield,
                                    None,
                                ),
                                target,
                            ),
                            make_component_variation(
                                bkg_const_hist,
                                np.clip(
                                    const_nominal * (1.0 - bkg_amp_step_fraction),
                                    floor_yield,
                                    None,
                                ),
                                target,
                            ),
                        ],
                        name=f"bkgC_eta{ieta}_pt{ipt}_charge{icharge}",
                        process="bkgConst",
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

    _add_aem_systematics(writer, variation_hist, signal_hist, args, labels, n_eta_bins)

    meta_data_dict = {}
    if args.signalNormPOI:
        meta_data_dict["signal_process"] = args.signalProcess

    writer.write(args.outfolder, outname, meta_data_dict=meta_data_dict)


if __name__ == "__main__":
    main()
