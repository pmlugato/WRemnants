"""Write the rabbit tensor for the CVH NanoAOD B+ -> J/psi K fit.

Reads the `fitbins` histogram from `scripts/histmakers/btojpsik_cvhnano.py` and
writes a single-channel tensor with exactly two processes:

    BuJpsiK      true B+ -> J/psi K+ from the inclusive simulation, fixed shape
    BuJpsiKBkg   flat template, seeded from the data residual

and, with `--scaleVariations`, one constrained nuisance per (A/e/M, eta group)
taken from the histmaker's `nominal_muonScaleSyst_responseWeights` tensor. Every
other free parameter of the fit belongs to a rabbit param model, so the tensor
holds templates, data and the scale variations and nothing else. That is the
opposite of `btojpsik_tensor.py`, whose per-bin normalisations and background
amplitudes are template-morph NOIs, and it is why this tensor is small.

The variation histogram is pushed through **the same** `reduce_fitbins` call as
the nominal, not a parallel implementation of it. That is the whole point of
design decision 10: the pT overflow fold, the mass window and the pT rebin cannot
drift between the nominal and its variations because there is one code path.

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

    parser.add_argument(
        "--scaleVariations",
        action="store_true",
        help="Add the A/e/M momentum-scale nuisances from the histmaker's "
        "nominal_muonScaleSyst_responseWeights tensor. Requires the histmaker to "
        "have been run with --scaleVariations.",
    )
    parser.add_argument(
        "--scaleVarHist",
        default="nominal_muonScaleSyst_responseWeights",
        help="Name of the variation histogram in the simulated dataset.",
    )
    parser.add_argument(
        "--scaleParams",
        nargs="+",
        default=["A", "e", "M"],
        help="Parameter names in the order the histmaker packed them along the "
        "'unc' axis. make_jpsi_crctn_unc_helper writes "
        "unc = eta_group * nparams + iparam, so this order is what turns an index "
        "back into a name and must match it.",
    )
    parser.add_argument(
        "--scaleVarScale",
        type=float,
        default=1.0,
        help="Inflate or shrink every variation about the nominal. >1 makes the "
        "nuisances easier to see in a first fit; the physical size is 1.0.",
    )
    parser.add_argument(
        "--scaleVarTrim",
        type=float,
        default=0.0,
        help="Set a variation bin back to nominal where |var - nom| is below this. "
        "0 disables. Only worth using if empty-cell noise destabilises the fit.",
    )
    parser.add_argument(
        "--skipScaleParams",
        nargs="*",
        default=[],
        help="Parameter names to leave out of the tensor, e.g. M, whose variation "
        "is a few keV and which no realistic yield constrains.",
    )
    return parser.parse_args()


def split_scale_variations(
    hsyst, nominal, params, scale=1.0, trim=0.0, skip=(), logger=None
):
    """Per-nuisance up/down template pairs from a reduced variation histogram.

    `hsyst` carries the nominal's axes plus ('unc', 'downUpVar'); `nominal` is the
    signal template the variations belong to, already on the fit binning and
    already scaled. Both must have come through the same reduction.

    The 'unc' index runs eta group fastest-outer, parameter fastest-inner --
    `unc = group * nparams + iparam` -- because that is how
    `make_jpsi_crctn_unc_helper` fills `var_mat`. Getting this backwards would
    silently relabel A as e, which no downstream check would catch, so the
    convention is asserted against the axis size rather than assumed.
    """
    names = list(hsyst.axes.name)
    for required in ("unc", "downUpVar"):
        if required not in names:
            raise RuntimeError(
                f"The variation histogram has no '{required}' axis: {names}"
            )
    nunc = hsyst.axes["unc"].size
    nparams = len(params)
    if nunc % nparams:
        raise RuntimeError(
            f"{nunc} uncertainty bins is not a multiple of {nparams} parameters "
            f"{params}; --scaleParams does not describe this histogram"
        )
    ngroups = nunc // nparams

    # WRemnants writes downUpVar index 0 for down and 1 for up.
    DOWN, UP = 0, 1

    ups, downs, out_names, groups = [], [], [], []
    for iparam, param in enumerate(params):
        if param in skip:
            continue
        for group in range(ngroups):
            iunc = group * nparams + iparam
            up = hsyst[{"unc": iunc, "downUpVar": UP}].copy()
            down = hsyst[{"unc": iunc, "downUpVar": DOWN}].copy()
            for h in (up, down):
                if list(h.axes.name) != list(nominal.axes.name):
                    raise RuntimeError(
                        f"Variation axes {list(h.axes.name)} do not match the "
                        f"nominal's {list(nominal.axes.name)}"
                    )
                if scale != 1.0:
                    h.values()[...] = nominal.values() + scale * (
                        h.values() - nominal.values()
                    )
                if trim > 0.0:
                    close = np.abs(h.values() - nominal.values()) < trim
                    h.values()[close] = nominal.values()[close]
            ups.append(up)
            downs.append(down)
            out_names.append(f"scale{param}{group}")
            groups.append(f"scale{param}")

    if logger is not None:
        for name, up, down in zip(out_names, ups, downs):
            rel = np.abs(0.5 * (up.values() - down.values())) / np.maximum(
                nominal.values(), 1e-12
            )
            logger.info(
                "  %-12s largest bin shift %.3e, template-integral shift %.3e",
                name,
                float(rel.max()),
                float(
                    abs(0.5 * (up.values().sum() - down.values().sum()))
                    / max(nominal.values().sum(), 1e-12)
                ),
            )
    return ups, downs, out_names, groups


def main():
    args = parse_args()
    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    data_raw, mc_raw, (data_key, mc_key) = cvh.load_fitbins(
        args.infile, args.datasetData, args.datasetMC
    )
    logger.info("Data: %s   simulation: %s", data_key, mc_key)

    syst_raw = None
    if args.scaleVariations:
        syst_raw = cvh.load_named_hist(args.infile, mc_key, args.scaleVarHist)
        logger.info(
            "Scale variations from %s: axes %s",
            args.scaleVarHist,
            list(syst_raw.axes.name),
        )

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
    # The same call, with the same kwargs, on a histogram that carries two extra
    # trailing axes. Nothing in reduce_fitbins is written per-axis-count, so the
    # variations get the nominal's treatment by construction.
    syst_reduced = (
        cvh.reduce_fitbins(syst_raw, **reduce_kwargs) if syst_raw is not None else None
    )

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

    scale_names = []
    if syst_reduced is not None:
        # Same category selection and the same global scale as the nominal, so
        # that var/nom is preserved exactly and the nuisance is a pure shape (and
        # small-normalisation) variation rather than a rescaling.
        syst_signal = cvh.select_category(syst_reduced, args.signalCategory)
        syst_signal.view(flow=False)["value"] *= signal_scale
        syst_signal.view(flow=False)["variance"] *= signal_scale**2

        logger.info("Scale nuisances:")
        ups, downs, scale_names, groups = split_scale_variations(
            syst_signal,
            scaled_signal,
            args.scaleParams,
            scale=args.scaleVarScale,
            trim=args.scaleVarTrim,
            skip=args.skipScaleParams,
            logger=logger,
        )
        for up, down, name, group in zip(ups, downs, scale_names, groups):
            writer.add_systematic(
                [up, down],
                name,
                args.signalProcess,
                args.channel,
                groups=[group],
                constrained=True,
            )
        logger.info(
            "Added %d scale nuisances in %d groups: %s",
            len(scale_names),
            len(set(groups)),
            sorted(set(groups)),
        )

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
        "scale_nuisances": scale_names,
        "scale_var_scale": args.scaleVarScale if args.scaleVariations else None,
        "scale_params": args.scaleParams if args.scaleVariations else [],
        "scale_params_skipped": list(args.skipScaleParams),
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
