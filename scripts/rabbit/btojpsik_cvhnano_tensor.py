"""Write the rabbit tensor for the CVH NanoAOD B+ -> J/psi K fit.

Reads the `fitbins` histogram from `scripts/histmakers/btojpsik_cvhnano.py` and
writes a single-channel tensor with exactly two processes:

    BuJpsiK      true B+ -> J/psi K+ from the inclusive simulation, fixed shape
    BuJpsiKBkg   flat template, seeded from the data residual

and no systematics at all. Every free parameter of the fit belongs to a rabbit
param model, so the tensor holds templates and data and nothing else. That is the
opposite of `btojpsik_tensor.py`, whose per-bin normalisations and background
amplitudes are template-morph NOIs, and it is why this tensor is small.

`btojpsik_tensor.py` is deliberately untouched: it produced the validated 2018
fits, and half of it is A/e/M variation machinery that has no input in this
NanoAOD (no per-leg gen kinematics, gap G2).

The fit is run with two composed param models, for example

    rabbit_fit.py <tensor> \\
      --paramModel AxisNormModel btojpsik_cvhnano BuJpsiK \\
            fit_bachelor_pt,fit_bachelor_eta,fit_bachelor_charge \\
      --paramModel AxisExpModel btojpsik_cvhnano BuJpsiKBkg jointCvh_mass \\
            fit_bachelor_pt,fit_bachelor_eta,fit_bachelor_charge \\
            fit_bachelor_eta

where the trailing argument groups the exponential slope on eta.

Example:

    python scripts/rabbit/btojpsik_cvhnano_tensor.py \\
        -i /ceph/submit/data/user/p/pmlugato/mz/histograms/btojpsik_cvhnano_full.hdf5 \\
        -o /ceph/submit/data/user/p/pmlugato/mz/tensors/ -p full
"""

import os

import numpy as np

from rabbit import tensorwriter
from wremnants.postprocessing import rabbit_cvhnano_helpers as cvh
from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_DATA,
    GEN_MC_CATEGORIES,
    GEN_SIGNAL,
)
from wremnants.utilities import common, parsing
from wums import logging


def parse_args():
    analysis_label = common.analysis_label(os.path.basename(__file__))
    parser, _ = parsing.common_parser(analysis_label)
    parser.description = __doc__

    parser.add_argument(
        "-i", "--infile", required=True, help="Histmaker output holding 'fitbins'."
    )
    parser.add_argument(
        "--channel",
        default="btojpsik_cvhnano",
        help="Channel name stored in the tensor. Channel and process names are "
        "specific to this channel on purpose: rabbit identifies processes by "
        "name across channels, so a process called 'signal' here and in the Z "
        "channel would become one process when merge_inputs.py combines them.",
    )
    parser.add_argument(
        "--datasetData",
        default="Charmonium",
        help="Data dataset key, or a substring of it.",
    )
    parser.add_argument(
        "--datasetMC",
        default="InclusiveBToJpsiX",
        help="Simulated dataset key, or a substring of it.",
    )
    parser.add_argument("--signalProcess", default="BuJpsiK")
    parser.add_argument("--bkgProcess", default="BuJpsiKBkg")
    parser.add_argument(
        "--signalCategory",
        type=int,
        default=GEN_SIGNAL,
        help="genCategory bin holding the signal. An index rather than a "
        "definition, so that tightening the truth match in the histmaker (Phase "
        "3a task 9.2 adds a momentum-closure requirement) needs no change here.",
    )
    parser.add_argument("--dataCategory", type=int, default=GEN_DATA)

    parser.add_argument(
        "--massLow",
        type=float,
        default=5.20,
        help="Low edge of the fit window. Must be an edge of the stored mass axis.",
    )
    parser.add_argument(
        "--massHigh",
        type=float,
        default=5.35,
        help="High edge of the fit window. The default [5.20, 5.35] at 10 bins is "
        "15 MeV per bin, which resolves the 32.5 MeV signal FWHM over 2.2 bins "
        "and leaves four sideband bins below the peak and three above. The old "
        "analysis' [5.20, 5.40] is also accepted and gives more upper sideband at "
        "1.6 bins per FWHM.",
    )
    parser.add_argument("--massBins", type=int, default=10)
    parser.add_argument(
        "--signalSigma",
        type=float,
        default=0.0138,
        help="Signal core width, used only to check the mass binning resolves the "
        "peak. Measured in Phase 2 from the narrow+broad+linear fit.",
    )
    parser.add_argument(
        "--peakHalfWidth",
        type=float,
        default=0.03,
        help="Half-width of the peak region, used for the seeding scale and for "
        "the mass-peak column of the occupancy table. Not a fit boundary.",
    )

    parser.add_argument(
        "--ptBins",
        type=int,
        default=3,
        help="Number of bachelor pT quantile bins. Phase 2 put the ceiling at 3 to "
        "4 at the 100-events-per-cell floor.",
    )
    parser.add_argument(
        "--ptEdges",
        type=float,
        nargs="+",
        default=None,
        help="Explicit pT edges, overriding --ptBins. Must be edges of the stored axis.",
    )
    parser.add_argument(
        "--etaRebin",
        type=int,
        default=1,
        help="Integer rebin factor for the 48-bin eta axis. 2, 3, 4, 6 and 8 divide it.",
    )
    parser.add_argument(
        "--collapse",
        type=str,
        nargs="*",
        default=[],
        choices=[cvh.PT_AXIS, cvh.ETA_AXIS, cvh.CHARGE_AXIS],
        help="Cell axes to sum over before writing. This is how the inclusive "
        "configuration and the coarse diagnostics are produced from the same "
        "input file.",
    )

    parser.add_argument(
        "--signalScale",
        type=float,
        default=None,
        help="Global scale applied to the signal template so the fitted "
        "normalisations start near 1. Defaults to the data/simulation ratio in "
        "the peak region. One factor for every cell, deliberately: a per-cell "
        "seed would write the quantity being measured into the starting point.",
    )
    parser.add_argument(
        "--seedHeadroom",
        type=float,
        default=0.8,
        help="Fraction of a cell's data the scaled signal template may occupy at "
        "the seed. A signal template cannot exceed the data, so this is a "
        "physical bound on the seed scale rather than a tuning knob; it is "
        "applied at a low percentile over cells so one nearly-empty cell cannot "
        "drive the scale to zero.",
    )
    parser.add_argument(
        "--pseudoData",
        default="none",
        choices=["none", "mcTotal"],
        help="'mcTotal' replaces the data histogram with the simulation's total "
        "over the three simulated categories, scaled to the data yield, so the "
        "true signal content is known and a background model's yield bias is "
        "measurable. The process list is unchanged.",
    )
    parser.add_argument(
        "--systematicType", default="log_normal", choices=["log_normal", "normal"]
    )
    return parser.parse_args()


def main():
    args = parse_args()
    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    data_raw, mc_raw, (data_key, mc_key) = cvh.load_fitbins(
        args.infile, args.datasetData, args.datasetMC
    )
    logger.info("Data: %s   simulation: %s", data_key, mc_key)

    # The pT quantile edges come from the simulated signal, so they have to be
    # derived before the reduction that consumes them. Derived on the folded,
    # windowed signal projection: quantiles of a distribution that still
    # includes the mass sidebands would be quantiles of the background.
    mc_folded = cvh.strip_flow(cvh.fold_overflow(mc_raw, cvh.PT_AXIS))
    mc_signal_windowed = cvh.select_category(
        cvh.select_mass_window(
            mc_folded, args.massLow, args.massHigh, args.massBins, args.signalSigma
        ),
        args.signalCategory,
    )
    if args.ptEdges:
        pt_edges = list(args.ptEdges)
        logger.info("Using explicit pT edges %s", pt_edges)
    else:
        pt_edges = cvh.quantile_edges(mc_signal_windowed, args.ptBins, cvh.PT_AXIS)

    reduce_kwargs = dict(
        mass_low=args.massLow,
        mass_high=args.massHigh,
        mass_bins=args.massBins,
        signal_sigma=args.signalSigma,
        pt_edges=pt_edges,
        eta_rebin=args.etaRebin,
        collapse_axes=args.collapse,
    )
    data_reduced = cvh.reduce_fitbins(data_raw, **reduce_kwargs)
    mc_reduced = cvh.reduce_fitbins(mc_raw, **reduce_kwargs)

    data_hist = cvh.select_category(data_reduced, args.dataCategory)
    signal_hist = cvh.select_category(mc_reduced, args.signalCategory)
    mc_total_hist = cvh.sum_categories(mc_reduced, GEN_MC_CATEGORIES)

    if list(data_hist.axes.name) != list(signal_hist.axes.name):
        raise RuntimeError(
            f"Data axes {list(data_hist.axes.name)} do not match signal axes "
            f"{list(signal_hist.axes.name)}"
        )
    for label, h in (("data", data_hist), ("signal", signal_hist)):
        if h.variances() is None:
            raise RuntimeError(f"The {label} template lost its variances in reduction")

    cell_axes = cvh.cell_axis_names(signal_hist)
    n_cells = (
        int(np.prod([signal_hist.axes[a].size for a in cell_axes])) if cell_axes else 1
    )
    n_mass = signal_hist.axes[cvh.MASS_AXIS].size
    logger.info(
        "Cell axes %s -> %d cells, %d mass bins, %d data bins",
        cell_axes or ["(none, inclusive)"],
        n_cells,
        n_mass,
        n_cells * n_mass,
    )

    # --- occupancy: the criterion the binning was chosen against ------------
    table = cvh.occupancy_table(
        mc_total_hist, signal_hist, cvh.peak_mask(signal_hist, args.peakHalfWidth)
    )
    print()
    print(cvh.format_occupancy(table))
    print()

    # --- seeding ------------------------------------------------------------
    # Sideband-subtracted, not a bare peak-region ratio: the data peak holds
    # signal plus background while the template holds signal alone, so an
    # unsubtracted ratio is too big by (1 + B/S) and starts most cells with a
    # signal template larger than their own data.
    seed_cap = None
    if args.signalScale is not None:
        signal_scale = args.signalScale
        logger.info(
            "Using the signal scale given on the command line: %.6g", signal_scale
        )
    else:
        signal_scale = cvh.signal_seed_scale(data_hist, signal_hist, args.peakHalfWidth)
        signal_scale, seed_cap = cvh.cap_seed_scale(
            signal_scale, data_hist, signal_hist, headroom=args.seedHeadroom
        )
    scaled_signal = signal_hist.copy()
    scaled_signal.view(flow=False)["value"] *= signal_scale
    scaled_signal.view(flow=False)["variance"] *= signal_scale**2

    if args.pseudoData == "mcTotal":
        # Substituting the observation, not adding processes: the fit model stays
        # the same two templates, and the true signal content is known, so a
        # background model's yield bias is measurable rather than inferred.
        data_yield = float(data_hist.values().sum())
        mc_yield = float(mc_total_hist.values().sum())
        pseudo_scale = data_yield / mc_yield if mc_yield > 0 else 1.0
        observed = mc_total_hist.copy()
        observed.view(flow=False)["value"] *= pseudo_scale
        observed.view(flow=False)["variance"] *= pseudo_scale**2
        true_signal = float(signal_hist.values().sum()) * pseudo_scale
        logger.info(
            "Pseudo-data from the simulated total, scaled by %.4f to the data "
            "yield %.6g. True signal content: %.6g events, i.e. a fitted "
            "normalisation of %.4f against this template.",
            pseudo_scale,
            data_yield,
            true_signal,
            true_signal / max(float(scaled_signal.values().sum()), 1e-12),
        )
    else:
        observed = data_hist
        pseudo_scale = None
        true_signal = None

    # --- background seed ----------------------------------------------------
    # Flat at the residual level per cell, so that AxisExpModel's default
    # lnAmpl = 0 starts the amplitude where the data is rather than at one event
    # per bin. The shape is then owned entirely by the param model.
    bkg_hist = signal_hist.copy()
    bkg_hist.view(flow=False)["value"] = 0.0
    bkg_hist.view(flow=False)["variance"] = 0.0

    mass_index = list(signal_hist.axes.name).index(cvh.MASS_AXIS)
    residual = observed.values() - scaled_signal.values()
    residual_per_cell = np.clip(residual.sum(axis=mass_index, keepdims=True), 0.0, None)
    # The floor is a small fraction of the cell's own data rather than an
    # absolute epsilon. A cell whose residual is zero still has data in it, so
    # starting its background at 1e-6 would put the optimiser eight orders of
    # magnitude from where it has to end up, in exactly the sparse cells that are
    # hardest to fit.
    observed_per_cell = observed.values().sum(axis=mass_index, keepdims=True)
    floor = np.maximum(0.05 * observed_per_cell, 1e-6)
    seed = np.maximum(residual_per_cell, floor) / n_mass
    bkg_hist.view(flow=False)["value"] = np.broadcast_to(seed, bkg_hist.values().shape)

    n_floored = int(np.sum(residual_per_cell <= 0.0))
    logger.info(
        "Background seed: %d of %d cells (%.1f%%) hit the positivity floor.",
        n_floored,
        n_cells,
        100.0 * n_floored / n_cells,
    )
    if n_floored > 0.1 * n_cells:
        # A handful of cells with no room for background is statistics. A large
        # fraction is the signal scale being too big, which is what an
        # unsubtracted peak-region ratio produces, and it starts the optimiser
        # somewhere it should not be.
        logger.warning(
            "%.0f%% of cells have data fully accounted for by the scaled signal. "
            "That is not purity, it is a signal scale that is too large: check "
            "the seed scale %.4f against the peak purity reported above.",
            100.0 * n_floored / n_cells,
            signal_scale,
        )

    # --- tensor -------------------------------------------------------------
    writer = tensorwriter.TensorWriter(systematic_type=args.systematicType)
    writer.add_channel(signal_hist.axes, name=args.channel)
    writer.add_data(observed, channel=args.channel)
    writer.add_process(scaled_signal, args.signalProcess, args.channel, signal=True)
    writer.add_process(bkg_hist, args.bkgProcess, args.channel)

    meta = {
        "signal_scale": signal_scale,
        "signal_scale_cap": seed_cap,
        "seed_cells_floored": n_floored,
        "signal_process": args.signalProcess,
        "bkg_process": args.bkgProcess,
        "signal_template_integral": float(scaled_signal.values().sum()),
        "unscaled_signal_integral": float(signal_hist.values().sum()),
        "mass_low": args.massLow,
        "mass_high": args.massHigh,
        "mass_bins": args.massBins,
        "signal_sigma": args.signalSigma,
        "pt_edges": [float(e) for e in pt_edges],
        "eta_rebin": args.etaRebin,
        "collapsed_axes": list(args.collapse),
        "cell_axes": cell_axes,
        "n_cells": n_cells,
        "data_dataset": data_key,
        "mc_dataset": mc_key,
        "pseudodata": args.pseudoData,
    }
    if true_signal is not None:
        meta["pseudodata_scale"] = pseudo_scale
        meta["pseudodata_true_signal"] = true_signal

    outname = os.path.splitext(os.path.basename(__file__))[0]
    if args.postfix:
        outname += f"_{args.postfix}"
    writer.write(args.outfolder, outname, meta_data_dict=meta)

    # The yield arithmetic, printed once so it is never reconstructed by hand.
    print()
    print("Recovering a yield from a fit result:")
    print(
        f"  signal template integral (scaled)   {meta['signal_template_integral']:.6g}"
    )
    print(f"  global signal scale applied         {signal_scale:.6f}")
    print(
        f"  simulated signal before scaling     {meta['unscaled_signal_integral']:.6g}"
    )
    print("  yield in a cell = fitted norm x that cell's scaled template integral")
    print("  total yield     = sum over cells; a fitted norm of 1 reproduces the")
    print("                    seed, which is the peak-region data/MC ratio")
    if true_signal is not None:
        print(f"  pseudo-data true signal             {true_signal:.6g}")
    print()
    logger.info("Wrote %s", os.path.join(args.outfolder, f"{outname}.hdf5"))


if __name__ == "__main__":
    main()
