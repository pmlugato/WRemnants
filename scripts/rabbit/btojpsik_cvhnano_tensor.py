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

import hist
import numpy as np

from rabbit import tensorwriter
from wremnants.postprocessing import rabbit_cvhnano_helpers as cvh
from wremnants.production.btojpsik_cvhnano_axes import (
    GEN_DATA,
    GEN_BD,
    GEN_BS,
    GEN_BU_KX,
    GEN_BU_PIX,
    GEN_MC_CATEGORIES,
    GEN_NON_SIGNAL_CATEGORIES,
    GEN_NO_B_ANCESTOR,
    GEN_OTHER_B,
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
        "--dataScale",
        type=float,
        default=1.0,
        help="Scale the observed histogram by this factor before writing the "
        "tensor. The SHAPE is untouched, so the model is wrong by exactly the "
        "same fraction, but every Poisson weight shrinks and a given fractional "
        "mismodelling costs sqrt(scale) fewer sigma. That isolates 'the model is "
        "wrong' from 'the model is wrong AND the data is precise enough to "
        "notice': if the fit converges at 0.01 and fails at 1.0 with the same "
        "shapes, the non-convergence is mis-specification amplified by "
        "statistics, not a broken parameterisation.",
    )
    parser.add_argument(
        "--resolutionVariation",
        type=float,
        default=0.0,
        help="Add a signal resolution nuisance by convolving the signal template "
        "along the mass axis with a Gaussian of this width in GeV; <=0 disables. "
        "This is the REFERENCE width: the nuisance is log-interpolated as "
        "(smeared/nominal)^theta, so the template only means 'a wider smear' "
        "when the fitted theta lands near 1. Use 0.005 for the per-eta primary "
        "(the data-minus-MC deficit runs 9 MeV at soft kaon pT to 3-4 MeV at "
        "hard pT, and an eta-only nuisance averages over that) and check the "
        "result at 0.003 (openspec add-cvhnano-per-eta-resolution, design 3).",
    )
    parser.add_argument(
        "--resolutionAxes",
        type=str,
        nargs="*",
        default=[],
        choices=[cvh.PT_AXIS, cvh.ETA_AXIS, cvh.CHARGE_AXIS],
        help="Cell axes over which the resolution nuisance varies. Empty (the "
        "default) keeps the single global 'sigRes'. 'fit_bachelor_eta' emits one "
        "free nuisance per eta bin of the FIT binning, each equal to the nominal "
        "outside its bin and to the smeared template inside it, named "
        "sigRes_eta{i}; 'fit_bachelor_eta fit_bachelor_pt' emits one per (eta, "
        "pT) group named sigRes_eta{i}_pt{j}. Eta-only is the physics choice: "
        "the multiple-scattering resolution is set by the material traversed, "
        "which is a function of eta, while the pT dependence is kinematic and "
        "already in the simulation. Requires --resolutionVariation > 0.",
    )
    parser.add_argument(
        "--bkgSeed",
        choices=["flat", "expSideband"],
        default="flat",
        help="Shape of the analytic background's TEMPLATE, i.e. where the optimiser "
        "starts. 'flat' (default, every recorded fit) spreads the per-cell "
        "residual evenly over the mass bins, so the exp*poly starts with zero "
        "slope against a spectrum that falls steeply: prefit loss 986k on the "
        "observation, and trust-krylov stalled in the indefinite-Hessian region "
        "it has to cross (2026-09-14). 'expSideband' fits exp(a + b*m) to each "
        "cell's sideband residual and uses that as the template. The fitted model "
        "class is unchanged -- a linear exponent is absorbed exactly by the "
        "poly's first-order term -- only the starting point moves.",
    )
    parser.add_argument(
        "--signalScaleReferenceYield",
        type=float,
        default=None,
        help="Metadata only: the data-fitted signal yield that --signalScale was "
        "derived from, so a pinned tensor records where its scale came from.",
    )
    parser.add_argument(
        "--resolutionConstrained",
        action="store_true",
        help="Give the resolution nuisance a unit Gaussian prior, so it behaves as "
        "a systematic and propagates into the A/e/M uncertainty. Default is "
        "UNCONSTRAINED: the MC width is measured to be wrong, so no prior width is "
        "defensible, and a constrained nuisance absorbing a known mismodelling is "
        "what produces multi-sigma pulls. Note the deficit is NOT symmetric -- the "
        "data falls off faster above the peak, and a low-side width ratio of 1.27 "
        "is worth chi2 761 on its own (openspec 11.3).",
    )
    parser.add_argument(
        "--mcBkgTemplate",
        choices=["none", "all", "comb", "split", "peak"],
        default="none",
        help="Add the simulation's own background shape as one or more template "
        "processes, alongside the analytic one. The analytic component is always "
        "kept, because the simulation has no prompt or non-b fake background and "
        "so cannot describe the data on its own.\n"
        "  none  -- analytic only (the original model)\n"
        "  all   -- one template from every non-signal MC category (1-6)\n"
        "  comb  -- one template from the no-b-ancestor (combinatorial) category 6\n"
        "  split -- two templates, combinatorial (6) and peaking-b (1-5), "
        "independently normalised, which tests whether the simulation got their "
        "ratio right\n"
        "  peak  -- the peaking-b template (1-5: B+ -> J/psi K X, B+ -> J/psi pi "
        "X, B0, Bs, other b) ALONE. The combinatorial MC template is 3.3x "
        "noisier than the data it would describe and degenerate with the "
        "analytic curvature (openspec aem 11.7, 9d.3), so the analytic component "
        "carries the smooth part and the template carries the B0 edge in the low "
        "shoulder that a polynomial cannot follow (11.6).",
    )
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
        "--massShiftB",
        type=int,
        default=0,
        help="Add an UNCONSTRAINED nuisance-of-interest 'massShiftB<N>MeV' whose "
        "up/down templates are the signal translated by +/- this many stored "
        "mass bins (5 MeV each). The B+ pole-mass parameter of interest, the "
        "analogue of the Z's massShiftZ100MeV. 0 = off.",
    )
    parser.add_argument(
        "--mcStatScale",
        type=float,
        default=1.0,
        help="Divide the variances of every simulation-derived template (signal, "
        "peaking MC, scale variations, resolution slabs) by this factor, as if "
        "that many times more events had been simulated. The template VALUES are "
        "untouched, so the Asimov dataset and the data-statistics term are "
        "unchanged; only the bin-by-bin MC-stat nuisances tighten. With a pinned "
        "--signalScale s the template variance per bin is s times the Asimov "
        "Poisson variance, so --mcStatScale s puts the simulation on the data's "
        "footing and 2s emulates twice as much simulation as data.",
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
        "--scaleVarParamScale",
        nargs="+",
        default=[],
        metavar="PARAM=FACTOR",
        help="Per-parameter stretch of the variations in log(var/nominal), "
        "applied after --scaleVarScale, e.g. 'A=5 e=5'. This changes the unit of the "
        "nuisance: build the histogram where the reweight is linear (1/1/1) and "
        "stretch here to the unit of the channel it is merged with "
        "(notes/aem_unit_mismatch.md 9.4). Recorded in the tensor metadata.",
    )
    parser.add_argument(
        "--scaleVarTrim",
        type=float,
        default=0.0,
        help="Set a variation bin back to nominal where |var - nom| is below this. "
        "0 disables. Only worth using if empty-cell noise destabilises the fit.",
    )
    parser.add_argument(
        "--scaleParamNaming",
        choices=["standalone", "joint"],
        default="standalone",
        help="How the A/e/M nuisances are named. 'standalone' (default) writes "
        "scale{param}{group} -- scaleA0, scalee0 -- which is what every existing "
        "B+ result was produced with. 'joint' writes {param}{group} -- A0, e0, M0 "
        "-- matching lowmass_tensor.py, so that merge_inputs.py correlates them "
        "with the Z and J/psi channels. merge_inputs correlates BY NAME AND "
        "NOTHING ELSE, so with the default naming a merged fit correlates zero "
        "parameters, converges, and looks entirely normal.",
    )
    parser.add_argument(
        "--jointEtaBins",
        type=int,
        default=None,
        help="Size of the shared nuisance eta basis, used ONLY to guard "
        "--scaleParamNaming joint. Must equal the number of eta groups the "
        "variation histogram actually carries, or the names would land on the "
        "wrong parameters. Required with --scaleParamNaming joint; the drivers "
        "pass it explicitly from their --nEta flag (no environment fallback).",
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
    hsyst,
    nominal,
    params,
    scale=1.0,
    trim=0.0,
    skip=(),
    logger=None,
    naming="standalone",
    joint_eta_bins=None,
    param_scale=None,
):
    """Per-nuisance up/down template pairs from a reduced variation histogram.

    ``param_scale`` maps a parameter name to a factor f applied in LOG space,
    up -> nominal * (up / nominal)**f, so that A and e can be stretched while M
    is not (see --scaleVarParamScale). rabbit interpolates in
    log(up / nominal), so this is the same nuisance with theta relabelled
    theta -> f * theta: identical slope at small theta, never negative.

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

    # The joint basis is shared by name, so a wrong count silently lands every
    # nuisance on a different eta bin -- and the fit still converges. Refuse.
    if naming == "joint":
        if joint_eta_bins is None:
            raise RuntimeError(
                "--scaleParamNaming joint needs --jointEtaBins"
            )
        if ngroups != joint_eta_bins:
            raise RuntimeError(
                f"--scaleParamNaming joint expects {joint_eta_bins} eta groups "
                f"but the variation histogram carries {ngroups}. Rerun the "
                f"histmaker with --scaleVariationEtaBins {joint_eta_bins}. "
                "Writing this tensor would correlate A/e/M with the wrong eta "
                "bins of the Z and J/psi channels, and nothing downstream would "
                "notice."
            )
        if skip:
            raise RuntimeError(
                f"--scaleParamNaming joint refuses --skipScaleParams {list(skip)}: "
                "a partially-correlated basis fits perfectly well and means "
                "something other than what it says."
            )

    # WRemnants writes downUpVar index 0 for down and 1 for up.
    DOWN, UP = 0, 1

    ups, downs, out_names, groups = [], [], [], []
    param_scale = param_scale or {}
    stretched = {}
    for iparam, param in enumerate(params):
        if param in skip:
            continue
        pfac = param_scale.get(param, 1.0)
        nkept = 0
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
                if pfac != 1.0:
                    nom, var = nominal.values(), h.values()
                    ok = (nom > 0) & (var > 0)
                    out = var.copy()
                    out[ok] = nom[ok] * (var[ok] / nom[ok]) ** pfac
                    # a bin with no positive nominal or variation carries no
                    # log-ratio to stretch; it is left as built
                    nkept += int(np.count_nonzero(~ok & (var != nom)))
                    h.values()[...] = out
                if trim > 0.0:
                    close = np.abs(h.values() - nominal.values()) < trim
                    h.values()[close] = nominal.values()[close]
            ups.append(up)
            downs.append(down)
            if naming == "joint":
                out_names.append(f"{param}{group}")
                groups.append(f"scale{param}")
            else:
                out_names.append(f"scale{param}{group}")
                groups.append(f"scale{param}")
        if pfac != 1.0:
            stretched[param] = (pfac, nkept)

    for param, (pfac, nkept) in stretched.items():
        if logger is not None:
            logger.info(
                "  %s variations stretched x%g in log(var/nominal); %d differing "
                "bins without a positive ratio left as built",
                param, pfac, nkept,
            )

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


def smear_along_mass(h, sigma, mass_axis_name, logger=None):
    """Convolve a histogram along its mass axis with a Gaussian of width `sigma`.

    Resolution is a convolution, so this is the variation's definition rather
    than an approximation of it. Three properties are enforced:

      * the integral is preserved per cell, so the nuisance is a pure SHAPE
        variation and cannot smuggle in a normalisation change (the signal
        normalisation is already a free parameter, and a variation that overlapped
        it would be degenerate with it);
      * the kernel is the exact bin-to-bin transfer for a density uniform within
        each source bin, so the width means GeV rather than bins **and the result
        is right when sigma is smaller than a bin**. The first version sampled
        the Gaussian at bin centres; at 20 MeV bins and a 5 MeV width that
        kernel is exp(-8) = 3e-4 off the diagonal, i.e. a delta function, and
        the "resolution nuisance" it produced moved the peak bin by 0.02%
        (found 2026-09-13 on the first per-eta tensor). The transfer matrix
        gives about 10% of each bin's content to each neighbour at that ratio,
        which is what a 5 MeV smear of a 20 MeV bin does;
      * edges are handled by renormalising each source bin's transfer over the
        bins that exist rather than by zero-padding, which would quietly leak
        yield out of the window.
    """
    from scipy.special import erf

    axes = list(h.axes.name)
    if mass_axis_name not in axes:
        raise ValueError(f"No '{mass_axis_name}' axis to smear along: {axes}")
    imass = axes.index(mass_axis_name)
    edges = np.asarray(h.axes[mass_axis_name].edges, dtype=float)
    width = float(edges[1] - edges[0])
    n = len(edges) - 1

    # T[i, j] = fraction of source bin j that lands in destination bin i, for a
    # density uniform over source bin j convolved with N(0, sigma). Closed form:
    # integrate Phi((E - x)/sigma) over x in the source bin, evaluated at both
    # destination edges. Uses G(t) = t*Phi(t) + phi(t), the antiderivative of
    # Phi, so that int_a^b Phi((E - x)/s) dx = s * [G((E - a)/s) - G((E - b)/s)].
    def G(t):
        Phi = 0.5 * (1.0 + erf(t / np.sqrt(2.0)))
        phi = np.exp(-0.5 * t * t) / np.sqrt(2.0 * np.pi)
        return t * Phi + phi

    a = edges[:-1][None, :]  # source bin low edges, shape (1, n)
    b = edges[1:][None, :]   # source bin high edges
    E = edges[:, None]       # destination edges, shape (n+1, 1)
    # C[k, j] = P(x + noise < E_k) integrated over source bin j, divided by width
    C = sigma * (G((E - a) / sigma) - G((E - b) / sigma)) / width
    T = C[1:, :] - C[:-1, :]  # (n, n): destination i, source j
    T = np.clip(T, 0.0, None)
    # Edge renormalisation: whatever would have left the window stays, spread
    # over the existing destination bins in proportion.
    T /= T.sum(axis=0, keepdims=True)

    out = h.copy()
    v = np.moveaxis(np.asarray(h.values(), dtype=float), imass, -1)
    var = np.moveaxis(np.asarray(h.variances(), dtype=float), imass, -1)
    sv = v @ T.T
    # Variances transform with the square of the kernel; this keeps the smeared
    # template's own MC statistical uncertainty honest rather than copying it.
    svar = var @ (T.T ** 2)
    out.view(flow=False)["value"] = np.moveaxis(sv, -1, imass)
    out.view(flow=False)["variance"] = np.moveaxis(svar, -1, imass)

    before = float(h.values().sum())
    after = float(out.values().sum())
    if before > 0 and abs(after / before - 1.0) > 1e-6 and logger is not None:
        logger.warning(
            "Smearing changed the integral by %.3g (%.6g -> %.6g); the kernel "
            "normalisation is meant to prevent this.",
            after / before - 1.0, before, after,
        )
    if logger is not None:
        logger.info(
            "Smeared the signal template by %.2f MeV (%.2f bins of %.2f MeV; an "
            "interior bin keeps %.1f%% of its content); integral %.6g -> %.6g "
            "(ratio %.9f)",
            sigma * 1e3, sigma / width, width * 1e3,
            100.0 * float(np.mean(np.diag(T)[1:-1])) if n > 2 else 100.0 * T[0, 0],
            before, after, after / max(before, 1e-12),
        )
    return out


def smear_at_stored_granularity(
    mc_raw, reduce_kwargs, category, signal_scale, sigma, nominal, logger=None,
    var_scale=None,
):
    """Smear the signal template at the histmaker's 5 MeV mass granularity, then
    rebin to the fit binning.

    The transfer matrix in `smear_along_mass` is exact for a density that is
    uniform within each source bin. At the fit binning (20 MeV) that assumption
    is wrong for a peak whose core is 14 MeV wide: the density inside the peak
    bin is far from flat, and treating it as flat makes the smeared content
    change linearly in sigma where a convolution of the true density changes
    as sigma^2. Measured 2026-09-14 on the first campaign: the fitted resolution
    uncertainty scaled as 5/3 between 5 and 3 MeV reference widths instead of
    (5/3)^2. The stored histogram has 5 MeV bins, so the smear is applied there
    -- where uniform-within-bin is a good approximation -- and the result is
    rebinned by the same factor the nominal was. Same reduction path
    (`reduce_fitbins`), same category selection, same scale, so the only
    difference from the nominal is the convolution.
    """
    stored = float(np.diff(mc_raw.axes[cvh.MASS_AXIS].edges)[0])
    low, high = reduce_kwargs["mass_low"], reduce_kwargs["mass_high"]
    n_fine = int(round((high - low) / stored))
    n_coarse = int(reduce_kwargs["mass_bins"])
    if n_fine % n_coarse:
        raise RuntimeError(
            f"{n_fine} stored mass bins in the window do not divide by "
            f"{n_coarse} fit bins; cannot smear at stored granularity"
        )
    factor = n_fine // n_coarse
    fine_kwargs = dict(reduce_kwargs, mass_bins=n_fine)
    fine = cvh.select_category(cvh.reduce_fitbins(mc_raw, **fine_kwargs), category)
    fine.view(flow=False)["value"] *= signal_scale
    fine.view(flow=False)["variance"] *= (
        signal_scale**2 if var_scale is None else var_scale
    )
    smeared_fine = smear_along_mass(fine, sigma, cvh.MASS_AXIS, logger=logger)
    smeared = smeared_fine[{cvh.MASS_AXIS: hist.rebin(factor)}]
    # The rebinned un-smeared fine template must be the nominal, exactly. If it
    # is not, the two took different paths and the variation is not a smear.
    check = fine[{cvh.MASS_AXIS: hist.rebin(factor)}]
    if not np.allclose(check.values(), nominal.values(), rtol=1e-9, atol=1e-9):
        raise RuntimeError(
            "Rebinned fine signal template differs from the nominal; the smear "
            "would not be a pure convolution of the template being fitted"
        )
    if logger is not None:
        logger.info(
            "Resolution smear applied at the stored %.1f MeV granularity (%d bins) "
            "and rebinned by %d to the %d fit bins",
            stored * 1e3, n_fine, factor, n_coarse,
        )
    return smeared


def shift_at_stored_granularity(
    mc_raw, reduce_kwargs, category, signal_scale, n_shift, nominal, logger=None,
    var_scale=None,
):
    """Translate the signal template along the mass axis by +/- n_shift stored
    bins (5 MeV each) and rebin to the fit binning: the B+ pole-mass variation.

    A pole-mass shift moves every reconstructed candidate by the same delta (no
    usable natural width, resolution 20-30 MeV), so the variation is a pure
    translation. It is done at the stored granularity, where it is exact, over
    a window n_shift stored bins wider than the fit window so that nothing is
    lost at the edges; the shifted fine template is then cropped to the fit
    window and rebinned by the same factor as the nominal. Returns (up, down),
    up = mass increased by n_shift x stored width.
    """
    stored = float(np.diff(mc_raw.axes[cvh.MASS_AXIS].edges)[0])
    low, high = reduce_kwargs["mass_low"], reduce_kwargs["mass_high"]
    n_fine = int(round((high - low) / stored))
    n_coarse = int(reduce_kwargs["mass_bins"])
    if n_fine % n_coarse:
        raise RuntimeError("stored bins do not divide the fit bins; cannot shift")
    factor = n_fine // n_coarse
    pad = int(n_shift)
    wide_kwargs = dict(
        reduce_kwargs,
        mass_low=low - pad * stored, mass_high=high + pad * stored,
        mass_bins=n_fine + 2 * pad, check_totals=False,
    )
    wide = cvh.select_category(cvh.reduce_fitbins(mc_raw, **wide_kwargs), category)
    wide.view(flow=False)["value"] *= signal_scale
    wide.view(flow=False)["variance"] *= (
        signal_scale**2 if var_scale is None else var_scale
    )
    imass = list(wide.axes.name).index(cvh.MASS_AXIS)
    v = np.asarray(wide.values(), dtype=float)
    var = np.asarray(wide.variances(), dtype=float)
    sl = [slice(None)] * v.ndim
    out = []
    for sign in (+1, -1):
        # rolling +pad moves content to higher mass; the wrapped edge bins fall
        # outside the fit window and are dropped by the crop below
        sv = np.roll(v, sign * pad, axis=imass); svar = np.roll(var, sign * pad, axis=imass)
        sl[imass] = slice(pad, pad + n_fine)
        fine = nominal.copy()  # coarse template, to be overwritten after rebinning
        fine_h = wide[{cvh.MASS_AXIS: slice(pad, pad + n_fine)}]
        fine_h.view(flow=False)["value"] = sv[tuple(sl)]
        fine_h.view(flow=False)["variance"] = svar[tuple(sl)]
        coarse = fine_h[{cvh.MASS_AXIS: hist.rebin(factor)}]
        fine.view(flow=False)["value"] = coarse.values()
        fine.view(flow=False)["variance"] = coarse.variances()
        out.append(fine)
    # closure: the unshifted crop must be the nominal
    sl[imass] = slice(pad, pad + n_fine)
    check = wide[{cvh.MASS_AXIS: slice(pad, pad + n_fine)}][{cvh.MASS_AXIS: hist.rebin(factor)}]
    if not np.allclose(check.values(), nominal.values(), rtol=1e-9, atol=1e-9):
        raise RuntimeError("cropped wide template differs from the nominal; the shift is not a pure translation")
    if logger is not None:
        up, down = out
        n0 = float(nominal.values().sum())
        logger.info(
            "B+ mass shift templates: +/-%d stored bins = +/-%.1f MeV, at %.1f MeV granularity, "
            "rebinned by %d; integral up/nominal %.6f, down/nominal %.6f (edge leakage only)",
            pad, pad * stored * 1e3, stored * 1e3, factor,
            float(up.values().sum()) / n0, float(down.values().sum()) / n0,
        )
    return out


def exp_sideband_seed(residual, cell_total, ref_hist, mass_index, peak_half_width):
    """Per-cell exponential template from the sideband residual.

    For each cell, weighted least squares of log(residual) against mass over
    the sideband bins with positive residual (weights = residual, the Poisson
    variance of a count), then exp(a + b (m - m0)) normalised to `cell_total`.
    Cells with fewer than three usable sideband bins keep the flat seed. Returns
    the shaped array, the number of fitted cells, the number of flat fallbacks,
    and the fitted slopes.
    """
    mass_axis = ref_hist.axes[cvh.MASS_AXIS]
    centres = np.asarray(mass_axis.centers, dtype=float)
    n_mass = len(centres)
    side = ~cvh.peak_mask(ref_hist, peak_half_width)
    m0 = float(centres.mean())
    x = centres - m0

    r = np.moveaxis(np.asarray(residual, dtype=float), mass_index, -1)
    tot = np.moveaxis(np.asarray(cell_total, dtype=float), mass_index, -1)
    shape = r.shape
    r2 = r.reshape(-1, n_mass)
    tot2 = tot.reshape(-1)
    out = np.empty_like(r2)
    n_fit = n_flat = 0
    slopes = []
    for i in range(r2.shape[0]):
        ok = side & (r2[i] > 0.0)
        if ok.sum() < 3:
            out[i] = tot2[i] / n_mass
            n_flat += 1
            continue
        w = r2[i][ok]
        A = np.column_stack([np.ones(ok.sum()), x[ok]])
        coef, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], np.log(w) * np.sqrt(w), rcond=None)
        shape_m = np.exp(coef[0] + coef[1] * x)
        out[i] = shape_m * (tot2[i] / shape_m.sum())
        slopes.append(float(coef[1]))
        n_fit += 1
    return np.moveaxis(out.reshape(shape), -1, mass_index), n_fit, n_flat, slopes


def split_resolution_by_axes(nominal, smeared, axes, logger=None):
    """Cut one smeared template into per-cell-group nuisances.

    Returns a list of (name, template) pairs, one per bin of the flattened index
    over `axes`. Template k equals `nominal` everywhere except in the cells of
    group k, where it equals `smeared`. Three properties are asserted on every
    slab, because a slab that fails any of them is a nuisance that does not mean
    what its name says:

      * its integral equals the nominal's (relative 1e-9): the convolution
        preserves the integral per cell, so the slab must too;
      * it is bit-identical to the nominal outside its group;
      * at least one bin inside its group moved.

    Names carry the bin indices, `sigRes_eta7` or `sigRes_eta7_pt1`, so a name
    can be read back to a cell without the tensor.
    """
    names_all = list(nominal.axes.name)
    for a in axes:
        if a not in names_all:
            raise ValueError(f"No '{a}' axis to split the resolution over: {names_all}")
    idx = [names_all.index(a) for a in axes]
    sizes = [nominal.axes[a].size for a in axes]
    short = {cvh.ETA_AXIS: "eta", cvh.PT_AXIS: "pt", cvh.CHARGE_AXIS: "q"}

    nom_v = np.asarray(nominal.values(), dtype=float)
    nom_var = np.asarray(nominal.variances(), dtype=float)
    sm_v = np.asarray(smeared.values(), dtype=float)
    sm_var = np.asarray(smeared.variances(), dtype=float)
    nom_total = float(nom_v.sum())

    out = []
    for flat in range(int(np.prod(sizes))):
        multi = np.unravel_index(flat, sizes)
        mask = np.ones(nom_v.shape, dtype=bool)
        for ax_i, b in zip(idx, multi):
            sel = [slice(None)] * nom_v.ndim
            sel[ax_i] = slice(b, b + 1)
            keep = np.zeros(nom_v.shape, dtype=bool)
            keep[tuple(sel)] = True
            mask &= keep
        h = nominal.copy()
        v = np.where(mask, sm_v, nom_v)
        var = np.where(mask, sm_var, nom_var)
        h.view(flow=False)["value"] = v
        h.view(flow=False)["variance"] = var

        name = "sigRes_" + "_".join(
            f"{short.get(a, a)}{b}" for a, b in zip(axes, multi)
        )
        total = float(v.sum())
        if nom_total > 0 and abs(total / nom_total - 1.0) > 1e-9:
            raise RuntimeError(
                f"{name}: integral {total:.9g} != nominal {nom_total:.9g}"
            )
        if np.any(v[~mask] != nom_v[~mask]):
            raise RuntimeError(f"{name}: differs from the nominal outside its slab")
        moved = np.abs(v[mask] - nom_v[mask])
        if moved.size == 0 or moved.max() <= 0.0:
            raise RuntimeError(
                f"{name}: no bin inside the slab moved; the slab is empty or the "
                "smear did nothing there"
            )
        if logger is not None:
            peak = float(nom_v[mask].max()) if mask.any() else 0.0
            logger.info(
                "  %-18s cells %5d  largest single-bin shift %.4g events (%.2f%% of "
                "the slab's peak bin)",
                name, int(mask.sum()), float(moved.max()),
                100.0 * float(moved.max()) / max(peak, 1e-12),
            )
        out.append((name, h))
    return out


def main():
    args = parse_args()
    param_scale = {}
    for item in args.scaleVarParamScale:
        name, _, factor = item.partition("=")
        if name not in args.scaleParams or not factor:
            raise ValueError(
                f"--scaleVarParamScale entry '{item}' must be PARAM=FACTOR with "
                f"PARAM in {args.scaleParams}"
            )
        param_scale[name] = float(factor)
    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    # Refused rather than ignored: a per-eta resolution with no reference width
    # is a request for nothing, and silently emitting nothing is the failure
    # mode the flags exist to remove.
    if args.resolutionAxes and args.resolutionVariation <= 0.0:
        raise RuntimeError(
            f"--resolutionAxes {args.resolutionAxes} needs --resolutionVariation > 0"
        )
    collapsed_res = [a for a in args.resolutionAxes if a in args.collapse]
    if collapsed_res:
        raise RuntimeError(
            f"--resolutionAxes names {collapsed_res}, which --collapse sums over; "
            "a nuisance cannot vary over an axis the tensor does not have"
        )

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
    if args.dataScale != 1.0:
        # Applied HERE, before the signal scale, the background seed and the MC
        # template scaling are derived, so every one of them is computed from the
        # scaled data and the fit starts at a consistent point. Applying it later
        # -- after the signal template has been scaled to the FULL data yield --
        # leaves the model over-predicting by 1/dataScale at iteration zero, and
        # then the test measures the optimiser's ability to travel a factor of
        # 100 rather than anything about statistics.
        data_hist = data_hist.copy()
        data_hist.view(flow=False)["value"] *= args.dataScale
        data_hist.view(flow=False)["variance"] *= args.dataScale
        logger.info(
            "Data scaled by %.4g -> %.6g events. Shape unchanged, so the model is "
            "wrong by the same FRACTION, but per-bin precision is %.2fx worse and "
            "a given mismodelling costs sqrt(%.4g) fewer sigma.",
            args.dataScale, float(data_hist.values().sum()),
            1.0 / np.sqrt(args.dataScale), args.dataScale,
        )
    signal_hist = cvh.select_category(mc_reduced, args.signalCategory)
    mc_total_hist = cvh.sum_categories(mc_reduced, GEN_MC_CATEGORIES)

    # MC background templates. Names are fixed rather than derived so the fit
    # script can name them on the --paramModel command line without knowing which
    # variant was built.
    PEAKING = [GEN_BU_KX, GEN_BU_PIX, GEN_BD, GEN_BS, GEN_OTHER_B]
    mc_bkg_specs = {
        "none": [],
        "all": [("BuJpsiKMCBkg", GEN_NON_SIGNAL_CATEGORIES)],
        "comb": [("BuJpsiKMCComb", [GEN_NO_B_ANCESTOR])],
        "split": [
            ("BuJpsiKMCComb", [GEN_NO_B_ANCESTOR]),
            ("BuJpsiKMCPeak", PEAKING),
        ],
        "peak": [("BuJpsiKMCPeak", PEAKING)],
    }[args.mcBkgTemplate]
    mc_bkg_hists = []
    for name, cats in mc_bkg_specs:
        h = cvh.sum_categories(mc_reduced, cats)
        mc_bkg_hists.append((name, cats, h))
        logger.info(
            "MC background template %s from categories %s: %.6g raw events",
            name,
            cats,
            float(h.values().sum()),
        )

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
    if args.mcStatScale <= 0.0:
        raise ValueError("--mcStatScale must be positive")
    # Every simulation-derived template gets its values scaled by signal_scale
    # and its variances by var_scale. With mcStatScale = 1 the two are the usual
    # pair (s, s^2); a larger mcStatScale shrinks the variances alone.
    var_scale = signal_scale**2 / args.mcStatScale
    if args.mcStatScale != 1.0:
        logger.info(
            "MC-stat scale %.4g: template variances divided by %.4g, values untouched. "
            "Template variance / Asimov Poisson variance per bin is %.3f (was %.3f).",
            args.mcStatScale, args.mcStatScale,
            signal_scale / args.mcStatScale, signal_scale,
        )
    scaled_signal = signal_hist.copy()
    scaled_signal.view(flow=False)["value"] *= signal_scale
    scaled_signal.view(flow=False)["variance"] *= var_scale

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
    # The MC templates carry the same global scale as the signal, so a fitted
    # normalisation of 1 means "the simulation's own prediction" and the fit is a
    # measurement of the composition rather than of an arbitrary scale. Their
    # variances are KEPT -- unlike the analytic component's, these are genuine MC
    # statistical uncertainties and deserve bin-by-bin nuisances.
    scaled_mc_bkgs = []
    for name, cats, h in mc_bkg_hists:
        sh = h.copy()
        sh.view(flow=False)["value"] *= signal_scale
        sh.view(flow=False)["variance"] *= var_scale
        scaled_mc_bkgs.append((name, sh))
        logger.info(
            "  %s scaled by %.4f -> %.6g events (%.2f%% of the data window)",
            name,
            signal_scale,
            float(sh.values().sum()),
            100.0 * float(sh.values().sum()) / max(float(observed.values().sum()), 1e-12),
        )

    bkg_hist = signal_hist.copy()
    bkg_hist.view(flow=False)["value"] = 0.0
    bkg_hist.view(flow=False)["variance"] = 0.0

    mass_index = list(signal_hist.axes.name).index(cvh.MASS_AXIS)
    # Subtract the MC templates as well as the signal: the analytic component is
    # there for what the simulation does NOT contain (prompt and non-b fakes), so
    # seeding it from data-minus-signal alone would double count everything the
    # templates already describe.
    residual = observed.values() - scaled_signal.values()
    for _, sh in scaled_mc_bkgs:
        residual = residual - sh.values()
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
    if args.bkgSeed == "expSideband":
        shaped, n_fit, n_flat, slopes = exp_sideband_seed(
            residual, seed * n_mass, signal_hist, mass_index, args.peakHalfWidth
        )
        bkg_hist.view(flow=False)["value"] = shaped
        logger.info(
            "Background seed 'expSideband': exp(a + b m) fitted to the sideband "
            "residual in %d of %d cells (%d fell back to flat); slope b median %.2f, "
            "range [%.2f, %.2f] per GeV. Cell integrals unchanged.",
            n_fit, n_cells, n_flat,
            float(np.median(slopes)) if slopes else float("nan"),
            float(np.min(slopes)) if slopes else float("nan"),
            float(np.max(slopes)) if slopes else float("nan"),
        )

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
    for name, sh in scaled_mc_bkgs:
        writer.add_process(sh, name, args.channel)
    writer.add_process(bkg_hist, args.bkgProcess, args.channel)

    scale_names = []
    if syst_reduced is not None:
        # Same category selection and the same global scale as the nominal, so
        # that var/nom is preserved exactly and the nuisance is a pure shape (and
        # small-normalisation) variation rather than a rescaling.
        syst_signal = cvh.select_category(syst_reduced, args.signalCategory)
        syst_signal.view(flow=False)["value"] *= signal_scale
        syst_signal.view(flow=False)["variance"] *= var_scale

        logger.info("Scale nuisances:")
        ups, downs, scale_names, groups = split_scale_variations(
            syst_signal,
            scaled_signal,
            args.scaleParams,
            scale=args.scaleVarScale,
            trim=args.scaleVarTrim,
            skip=args.skipScaleParams,
            naming=args.scaleParamNaming,
            joint_eta_bins=args.jointEtaBins,
            logger=logger,
            param_scale=param_scale,
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

    # --- resolution nuisance -------------------------------------------------
    # ONE smeared template: the signal convolved with a Gaussian. rabbit's
    # mirror=True then interpolates multiplicatively as (up/nominal)^theta, so
    # theta = 0 is the nominal, +1 the smeared template, and -1 a narrower one
    # obtained by log-extrapolation. Positive for every theta, which a hand-built
    # additive mirror is not: 2*nominal - up went negative in 274 of 1280 bins,
    # because smearing by 1.78 bins moves enough yield into the tails that twice
    # the nominal cannot cover it there.
    #
    # A log-symmetric variation is also the right shape for this nuisance. The
    # data/template width ratio is multiplicative (1.168), not additive.
    #
    # With --resolutionAxes the one smeared template is cut into slabs: nuisance
    # i is the nominal everywhere except in cell-group i, where it is the smeared
    # template. Each slab is a pure shape variation (the convolution preserves
    # the integral per cell) and has no lever outside its group, so 24 of them
    # on eta are 24 independent widths that share nothing with the per-eta
    # signal normalisation.
    resolution_names = []
    if args.resolutionVariation > 0.0:
        res_up = smear_at_stored_granularity(
            mc_raw, reduce_kwargs, args.signalCategory, signal_scale,
            args.resolutionVariation, scaled_signal, logger=logger,
            var_scale=var_scale,
        )
        constrained_label = (
            "CONSTRAINED (unit Gaussian prior)" if args.resolutionConstrained
            else "UNCONSTRAINED (no prior; the data determines it)"
        )
        if not args.resolutionAxes:
            shift = np.abs(res_up.values() - scaled_signal.values()).max()
            logger.info(
                "Resolution nuisance 'sigRes': %s, largest single-bin shift %.4g "
                "events (%.2f%% of the peak bin)",
                constrained_label,
                shift, 100.0 * shift / max(float(scaled_signal.values().max()), 1e-12),
            )
            writer.add_systematic(
                res_up,
                "sigRes",
                args.signalProcess,
                args.channel,
                mirror=True,
                groups=["resolution"],
                constrained=bool(args.resolutionConstrained),
            )
            resolution_names.append("sigRes")
        else:
            slab_templates = split_resolution_by_axes(
                scaled_signal, res_up, args.resolutionAxes, logger=logger
            )
            logger.info(
                "Resolution nuisances over %s: %d slabs, %s",
                args.resolutionAxes, len(slab_templates), constrained_label,
            )
            for name, h in slab_templates:
                writer.add_systematic(
                    h,
                    name,
                    args.signalProcess,
                    args.channel,
                    mirror=True,
                    groups=["resolution"],
                    constrained=bool(args.resolutionConstrained),
                )
                resolution_names.append(name)

    mass_shift_name = None
    if args.massShiftB > 0:
        stored_mev = float(np.diff(mc_raw.axes[cvh.MASS_AXIS].edges)[0]) * 1e3
        mass_shift_name = f"massShiftB{int(round(args.massShiftB * stored_mev))}MeV"
        up, down = shift_at_stored_granularity(
            mc_raw, reduce_kwargs, args.signalCategory, signal_scale,
            args.massShiftB, scaled_signal, logger=logger, var_scale=var_scale,
        )
        writer.add_systematic(
            [up, down],
            mass_shift_name,
            args.signalProcess,
            args.channel,
            mirror=False,
            groups=["massShiftB"],
            constrained=False,
            noi=True,   # a nuisance-of-interest: rabbit computes impacts ON it, like massShiftZ100MeV
        )
        logger.info("Added the B+ mass POI '%s' (unconstrained, signal template only)", mass_shift_name)

    meta = {
        "signal_scale": signal_scale,
        "mass_shift_noi": mass_shift_name,
        "data_scale": args.dataScale,
        "mc_stat_scale": args.mcStatScale,
        "resolution_variation": args.resolutionVariation,
        "resolution_constrained": bool(args.resolutionConstrained),
        "bkg_seed": args.bkgSeed,
        "resolution_axes": list(args.resolutionAxes),
        "resolution_nuisances": resolution_names,
        "signal_scale_reference_yield": (
            args.signalScaleReferenceYield if args.signalScale is not None else None
        ),
        "signal_scale_cap": seed_cap,
        "seed_cells_floored": n_floored,
        "signal_process": args.signalProcess,
        "bkg_process": args.bkgProcess,
        "mc_bkg_template": args.mcBkgTemplate,
        "mc_bkg_processes": [n for n, _ in scaled_mc_bkgs],
        "mc_bkg_integrals": {n: float(h.values().sum()) for n, h in scaled_mc_bkgs},
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
        # Per-parameter factor on top of scale_var_scale. One unit of A_i in this
        # tensor is (histmaker width) x scale_var_scale x this; the joint-fit
        # guard (check_joint_tensor.py --bTensor) reads both.
        "scale_var_param_scale": param_scale if args.scaleVariations else {},
        "input_histogram": os.path.realpath(args.infile),
        "scale_params": args.scaleParams if args.scaleVariations else [],
        "scale_params_skipped": list(args.skipScaleParams),
        "scale_param_naming": args.scaleParamNaming,
        "joint_eta_bins": args.jointEtaBins if args.scaleParamNaming == "joint" else None,
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
