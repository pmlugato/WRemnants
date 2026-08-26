"""B+ -> J/psi K histmaker for the CVH NanoAOD.

The input is a NanoAOD produced from our TkAlJpsiX AlCaReco by
``Analysis/HitAnalyzer/test/runCvhBplusJpsiK.py``, carrying flat ``BuJpsiK_*``
candidates with three CVH fit arms. It shares no physics branch name with the
BMM-Tools NanoAOD that ``btojpsik.py`` reads, which is left untouched so the
validated 2018 results stay reproducible.

Scope is kinematic distributions and data/MC comparison. The kaon scale
variations are not wired: they need gen kaon pT, which this NanoAOD does not
provide (gap G2 in the proposal).

Example:

    python scripts/histmakers/btojpsik_cvhnano.py \
        --era 2016PostVFP --dataPath /ceph/submit/data/group/cms/store \
        --filterProcs Charmonium_2016PostVFP_cvhNano InclusiveBToJpsiX_2016PostVFP \
        --allaxes --cutflow --noScaleToData --theoryCorr none \
        -o /ceph/submit/data/user/p/pmlugato/mz/histograms/ -p v1
"""

import json
import os

import hist

import narf
from wremnants.production import btojpsik_cvhnano_selections as sel
from wremnants.production import muon_calibration
from wremnants.production.btojpsik_cvhnano_axes import (
    ARMS,
    CAND,
    DEFAULT_ARM,
    EVENT_SCALAR_AXES,
    GEN_SIGNAL,
    ID_AXES,
    MC_ONLY_AXES,
    all_axes,
    axis_ancestor_state,
    axis_fit_bachelor_charge,
    axis_fit_bachelor_eta,
    axis_fit_bachelor_pt,
    axis_gen_category,
    axis_gen_truth_eta,
    axis_gen_truth_pt,
    axis_id_bachelor_pt,
    axis_n_legs_matched,
    axis_rel_dpt,
    dimuon_vtx_arm,
)
from wremnants.production.datasets.dataset_tools import getDatasets
from wremnants.production.histmaker_tools import (
    aggregate_groups,
    scale_to_data,
    write_analysis_output,
)
from wremnants.utilities import common, parsing
from wums import logging

analysis_label = common.analysis_label(os.path.basename(__file__))
parser, initargs = parsing.common_parser(analysis_label)

parser.add_argument(
    "--arm",
    type=str,
    default=DEFAULT_ARM,
    choices=list(ARMS.keys()),
    help="Fit arm used for the nominal mass, the selections and candidate "
    "ranking. All arms are histogrammed regardless of this choice.",
)
parser.add_argument(
    "--trigger",
    type=str,
    default="DoubleMu4_JpsiTrk_Displaced",
    help="HLT path, applied identically to data and MC. The 2016 default of "
    "btojpsik.py (DoubleMu4_PsiPrimeTrk_Displaced) fires in 0.08%% of this "
    "data and must not be used. The low-threshold paths were prescaled in "
    "2016H but not in MC, so they are not comparable either.",
)
parser.add_argument("--allaxes", action="store_true", help="Histogram every known axis")
parser.add_argument("--axes", type=str, nargs="*", default=[], help="Axes to histogram")
parser.add_argument(
    "--selectionHists",
    action="store_true",
    help="Store the nominal histograms after each selection as well as at the end",
)
parser.add_argument("--cutflow", action="store_true", help="Print the cutflow table")
parser.add_argument(
    "--fitBinsHist",
    action="store_true",
    help="Book the 5D differential histogram (bachelor pT, eta, charge, mass, "
    "genCategory) used to measure whether the simulation can support the "
    "intended fit binning.",
)
parser.add_argument(
    "--qualityHists",
    action="store_true",
    help="Book one 3D histogram per bachelor quality variable: the variable "
    "against bachelor pT against genCategory. Needed for the kaon-ID study, "
    "because every candidate discriminant depends on pT -- displacement "
    "significance grows with it and dE/dx separates K from pi only below "
    "~0.7 GeV/c -- so an inclusive 1D distribution answers the wrong question.",
)
parser.add_argument(
    "--saveCutflow", type=str, default=None, help="Directory for the cutflow table"
)
parser.add_argument(
    "--skipFileCheck",
    action="store_true",
    help="Do not drop input files whose production JSON sidecar reports a failure",
)
parser.add_argument(
    "--skipTriggerCheck",
    action="store_true",
    help="Do not verify the trigger column is present in every input file. The "
    "HLT menu varies between files, so skipping this risks a mid-loop failure.",
)

# Quality gate. jointCvhOk == 1 is a convergence flag, not a quality verdict:
# it admits fits with chi-square up to 8e29 and mass error up to 4e5. Defaults
# are the payload-side convergence threshold for edmRef, and roughly an order
# of magnitude beyond the 95th percentile of the surviving population for the
# rest, so they cut the pathological tail without sculpting the bulk.
parser.add_argument("--gateEdmRef", type=float, default=1e-2, help="Max edmRef")
parser.add_argument("--gateChisq", type=float, default=1e6, help="Max fit chi-square")
parser.add_argument("--gateNdof", type=float, default=200.0, help="Max fit ndof")
parser.add_argument(
    "--gateMassErr", type=float, default=1.0, help="Max fitted mass err"
)

# Kinematic and geometric selections.
parser.add_argument("--muonMinPt", type=float, default=3.0)
parser.add_argument("--muonMaxEta", type=float, default=2.4)
parser.add_argument("--requireGlobalMuons", action="store_true", default=False)
parser.add_argument(
    "--bachelorMinPt",
    type=float,
    default=1.0,
    help="Bachelor pT floor, default 1 GeV. Justified twice over. It makes the "
    "combinatorial describable by a monotonic analytic shape (exponential "
    "chi2/ndf 745.8 -> 41.2, yield bias 2.8x -> 1.11x, 96.8% of the signal "
    "retained), and it keeps the mass-aware e variation in the regime where "
    "the exact finite shift is well behaved -- the margin E - m_K against the "
    "e uncertainty is 22 MeV at pT 0.15 GeV and 622 MeV at 1 GeV. It also "
    "matches the production-level cut the 2018 sample already had, which makes "
    "the two comparable. Lower it to reach the soft region: the AlCaReco floor "
    "is 0.1 and the observed minimum in the NanoAOD is 0.151 GeV.",
)
parser.add_argument("--bachelorMaxEta", type=float, default=2.4)
parser.add_argument(
    "--bachelorMaxPt",
    type=float,
    default=-1.0,
    help="Bachelor pT ceiling; <0 disables. The 2018 analysis capped the kaon at "
    "8 GeV.",
)
parser.add_argument(
    "--scaleVariations",
    action="store_true",
    help="Book the A/e/M momentum-scale variations on the simulated template.",
)
parser.add_argument(
    "--scaleVarLegs",
    choices=["all", "bachelor"],
    default="all",
    help="Which legs the A/e/M variations are applied to. 'all' (default) "
    "varies bachelor, mu0 and mu1, because the muon A/e/M are real parameters "
    "of the same calibration and omitting them would understate the "
    "candidate's total scale uncertainty. 'bachelor' isolates the bachelor's "
    "contribution, so the muon contribution can be measured rather than "
    "assumed small.",
)
parser.add_argument(
    "--eConvention",
    choices=["pt", "energy"],
    default="pt",
    help="How the correction file's e parameter is interpreted. 'pt' (default) "
    "treats it as a shift of transverse momentum, which is how it was fitted: "
    "the central-value corrector applies (1 + A - e*k)*k with no cosh(eta), so "
    "any eta-dependent factor is already inside the fitted e_i(eta). 'energy' "
    "treats it as a shift of total energy, which carries an extra 1/cosh(eta) "
    "-- 1.0 at eta 0, 2.15 at 1.4, 5.56 at 2.4. See notes/E_IN_MATERIAL_TERM.md.",
)
parser.add_argument(
    "--scaleVarMuonSource",
    type=int,
    default=443,
    help="The network's per-leg class tag. 443 is the J/psi-leg class and is "
    "the default for every leg including the bachelor, following d0_mass.py -- "
    "there is no hadron class in the training. 1 is the prompt W/Z class; "
    "running with it measures how much the choice matters rather than "
    "assuming it does not.",
)
parser.add_argument(
    "--scaleVariationEtaBins",
    type=int,
    default=4,
    help="Number of eta bins the A/e/M variations are grouped into. The "
    "correction file has 48, which would give 144 scale nuisances on a channel "
    "with of order 1e4 signal candidates, so the default is deliberately "
    "coarse. Each fine eta bin inside a group is still moved by its own "
    "uncertainty under one shared nuisance -- a fully correlated variation "
    "within the group, not an average over it. Must divide 48. Set to 48 to "
    "recover the basis the joint mZ + J/psi fit correlates by name.",
)
parser.add_argument(
    "--jpsiFixedAUnc",
    type=float,
    default=None,
    help="Override the A uncertainty in every eta bin, as btojpsik.py does.",
)
parser.add_argument(
    "--jpsiFixedEUnc",
    type=float,
    default=None,
    help="Override the e uncertainty in every eta bin [GeV].",
)
parser.add_argument(
    "--jpsiFixedMUnc",
    type=float,
    default=None,
    help="Override the M uncertainty in every eta bin.",
)
parser.add_argument(
    "--genPtReweightSaturation",
    type=float,
    default=2.0,
    help="Gen-pt floor for the response network, in GeV; <=0 disables. The "
    "training sample floors gen pt at 2.0 GeV, so a bachelor below that is "
    "outside the model's domain even after the 1 GeV selection.",
)
parser.add_argument(
    "--genPtReweightSaturationMode",
    choices=["condition", "rescale"],
    default="condition",
    help="'condition' (default) floors only the network's conditioning input "
    "log pt_gen, leaving the response residual and the physical shift at the "
    "true kinematics. 'rescale' would scale both reco and gen pt, which "
    "preserves qopr but evaluates the shift at the saturated pt; it is not "
    "implemented here and is rejected by the argument check below.",
)
parser.add_argument("--massLow", type=float, default=5.0)
parser.add_argument("--massHigh", type=float, default=5.5)
parser.add_argument(
    "--minVtxProb",
    type=float,
    default=0.0,
    help="Candidate vertex probability floor. For the joint arm this is the "
    "chi-square upper-tail probability at the stored ndof, since that arm stores "
    "no vertex probability. The 2018 analysis cut this at 0.3, and it is the "
    "tightest cut in that selection.",
)
parser.add_argument("--minDimuonVtxProb", type=float, default=0.0)
parser.add_argument(
    "--minDimuonSl3d",
    type=float,
    default=0.0,
    help="Dimuon 3D decay-length significance. The 2018 cut was on the DIMUON's, "
    "not the candidate's; --minSl3d is the candidate's.",
)
parser.add_argument(
    "--maxDimuonAlphaBS",
    type=float,
    default=-1.0,
    help="Dimuon pointing angle to the beamspot; <0 disables. Again the 2018 cut "
    "was the dimuon's, not the candidate's --maxAlphaBS.",
)
parser.add_argument(
    "--minDimuonPt",
    type=float,
    default=0.0,
    help="Dimuon pT floor, on the dimuon rebuilt from the raw muon tracks. The "
    "2018 analysis cut at 7 GeV.",
)
parser.add_argument("--minSl3d", type=float, default=0.0)
parser.add_argument("--maxAlphaBS", type=float, default=-1.0, help="<0 disables")
parser.add_argument("--minLegsRefit", type=int, default=3)

# Bachelor and joint-fit quality. All default to off: each is under study in
# add-cvhnano-bachelor-quality-cuts and none is validated.
parser.add_argument(
    "--bachelorMinValidHits",
    type=int,
    default=-1,
    help="Bachelor valid-hit floor (strict >); <0 disables. Measured gain peaks "
    "near >13, not >10.",
)
parser.add_argument(
    "--minBachelorAbsD0",
    type=float,
    default=-1.0,
    help="Bachelor |d0| floor, PV-referenced; <0 disables.",
)
parser.add_argument(
    "--minBachelorAbsD0Sig",
    type=float,
    default=-1.0,
    help="Bachelor |d0|/sigma(d0) floor w.r.t. the primary vertex; <0 disables. "
    "The kaon ID: the background is a soft PROMPT track from the PV and the "
    "signal kaon is displaced. Implied peak purity 53%% at 1.0, 71%% at 1.5 and "
    "83%% at 2.0, for signal efficiencies 0.86, 0.80 and 0.74 -- against 58%% for "
    "kaon pT > 1 GeV, which also throws away the soft kaons.",
)
parser.add_argument(
    "--minBachelorDedx",
    type=float,
    default=-1.0,
    help="Bachelor harmonic-2 dE/dx floor; <0 disables. Mostly removes candidates "
    "whose estimator is unset rather than identifying kaons.",
)
parser.add_argument(
    "--maxJointCvhNiter",
    type=int,
    default=-1,
    help="Joint-fit iteration ceiling (inclusive); <0 disables. Strongest single "
    "discriminant measured, and the most dangerous: it is a diagnostic of the fit "
    "that produces the mass, so it can sculpt the peak.",
)
parser.add_argument(
    "--maxChisqNdof",
    type=float,
    default=-1.0,
    help="Joint-fit chi-square per degree of freedom ceiling; <0 disables. The "
    "normalised variable: ndof spans 32-66, so a raw chi-square cut varies in "
    "tightness with hit count and therefore with pT and eta.",
)

parser.add_argument(
    "--pileupWeights",
    type=str,
    default=None,
    help="JSON file of int(Pileup_nTrueInt) -> weight, applied to simulation "
    "only. Derive with scripts/plotting/btojpsik_cvhnano_pileup.py. Was keyed "
    "on the reconstructed vertex count while the NanoAOD carried no truth "
    "pileup (gap G4, now closed), which made the nPV agreement true by "
    "construction rather than a test.",
)
parser.add_argument(
    "--noForcedDecayWeight",
    action="store_true",
    help="Do not undo the generator's forced b-hadron decays. The forced BRs "
    "span 0.00047 to 0.911 across species, so leaving them in place gives a "
    "stacked simulation whose component fractions are not physical.",
)

parser.add_argument(
    "--triggerCheckCache",
    type=str,
    default=os.path.join(
        os.path.expanduser("~"), ".cache", "wremnants", "cvhnano_trigger_check.json"
    ),
    help="Where to remember that a given file list has already been checked for "
    "the trigger column. The check opens every input file and its answer does "
    "not change for a fixed production, so re-deriving it on every run is pure "
    "wall clock: it cost over an hour on 22 327 files. Set to an empty string to "
    "disable the cache.",
)
parser.add_argument(
    "--genDiagnostics",
    action="store_true",
    help="Book the generator-truth diagnostics: the gen-level B+ -> J/psi K+ "
    "denominator found by walking the Gen table with no reco requirement, the "
    "per-candidate matcher outcome (how many legs the producer accepted), and "
    "the bachelor pT response. Together these say whether the signal template "
    "is defined by a cut on the quantity being measured -- the per-leg match "
    "requires |dpt|/pt < 0.1, and the template is the population that passed "
    "it. Simulation only.",
)

parser = parsing.set_parser_default(parser, "aggregateGroups", [])
parser = parsing.set_parser_default(parser, "theoryCorr", [])
# The inclusive MC carries a placeholder cross section, so luminosity scaling
# would silently produce a wrong normalisation. Normalisation belongs to the
# tensor stage; default to leaving the yields raw.
parser = parsing.set_parser_default(parser, "noScaleToData", True)

args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

arm = ARMS[args.arm]
logger.info("Nominal fit arm: %s (%s)", arm.name, arm.label)

datasets = getDatasets(
    maxFiles=args.maxFiles,
    filt=args.filterProcs,
    excl=args.excludeProcs,
    base_path=args.dataPath,
    era=args.era,
    data_tags=[""],
    mc_tags=[""],
)

if not datasets:
    raise RuntimeError(
        "No datasets selected. The CVH NanoAOD samples are auxiliary, so they "
        "must be named explicitly with --filterProcs, e.g. "
        "--filterProcs Charmonium_2016PostVFP_cvhNano InclusiveBToJpsiX_2016PostVFP"
    )

# Drop truncated stubs left by killed production jobs before anything else
# touches the files: a file that exists is not a file that is complete.
if not args.skipFileCheck:
    for dataset in datasets:
        dataset.filepaths = sel.select_good_files(dataset.filepaths)
    datasets = [d for d in datasets if d.filepaths]

trigger_column = (
    args.trigger if args.trigger.startswith("HLT_") else f"HLT_{args.trigger}"
)
if not args.skipTriggerCheck:
    for dataset in datasets:
        sel.validate_trigger_available(
            dataset.filepaths, trigger_column, cache=args.triggerCheckCache
        )

# ---------------------------------------------------------------------------
# Momentum-scale variation helper
# ---------------------------------------------------------------------------

scale_unc_helper = None
if args.scaleVariations:
    if args.muonScaleVariation != "onnxReweight":
        raise ValueError(
            f"--muonScaleVariation {args.muonScaleVariation} is not wired in "
            "this channel. The analytic splines backend needs a per-leg "
            "response-weight column from a tflite response model, and the only "
            "kaon-trained one on disk predates this selection and this "
            "NanoAOD, so comparing against it now would confound three "
            "changes at once. It is a planned cross-check once that model is "
            "regenerated."
        )
    if args.genPtReweightSaturationMode != "condition":
        raise ValueError(
            "Only --genPtReweightSaturationMode condition is implemented here. "
            "'rescale' scales both reco and gen pt, which floors gen pt while "
            "preserving qopr but then evaluates the physical shift at the "
            "saturated pt; 'condition' keeps the shift honest and moves only "
            "the network's conditioning input."
        )
    if 48 % args.scaleVariationEtaBins:
        raise ValueError(
            f"--scaleVariationEtaBins {args.scaleVariationEtaBins} must divide "
            "the 48 eta bins of the correction file"
        )
    _, _, _, scale_unc_helper = muon_calibration.make_jpsi_crctn_helpers(
        common.calib_filepaths,
        muon_corr_mc=args.muonCorrMC,
        muon_corr_data=args.muonCorrData,
        scale_var_method=args.muonScaleVariation,
        scale_A=args.scale_A,
        scale_e=args.scale_e,
        scale_M=args.scale_M,
        make_uncertainty_helper=True,
        # No eigen-decomposition: the coarsened eta binning below is only
        # available without the covariance, and this channel constrains one
        # parameter group at a time rather than the full 48-bin basis.
        include_covariance=False,
        fixed_A_unc=args.jpsiFixedAUnc,
        fixed_e_unc=args.jpsiFixedEUnc,
        fixed_M_unc=args.jpsiFixedMUnc,
        variation_eta_bins=args.scaleVariationEtaBins,
        # Bachelor first, then the two muons -- SCALE_VAR_LEGS order. Explicit
        # per-leg rather than a one-element list, which both helpers broadcast
        # to every leg and would hand the muons a 494 MeV energy-loss term.
        reweight_mass=(
            [sel.KAON_MASS_GEV, 0.0, 0.0]
            if args.scaleVarLegs == "all"
            else [sel.KAON_MASS_GEV]
        ),
        # This channel applies no resolution smearing, so the model has to be
        # the one trained on un-smeared reco pt. Tied to --noSmearing rather
        # than chosen independently, so the two cannot drift apart.
        smearing=not args.noSmearing,
        cond_pt_gen_min=(
            args.genPtReweightSaturation if args.genPtReweightSaturation > 0 else None
        ),
        e_convention=args.eConvention,
    )
    if scale_unc_helper is None:
        raise ValueError(
            "--scaleVariations needs --muonCorrData massfit or lbl_massfit, "
            "since the uncertainties are read from that correction file"
        )
    logger.info(
        "Scale variations: %s, %d eta groups, %d nuisances, legs %s",
        args.muonScaleVariation,
        args.scaleVariationEtaBins,
        scale_unc_helper.tensor_axes[0].size,
        args.scaleVarLegs,
    )

for a in args.axes:
    if a not in all_axes:
        raise ValueError(
            f"{a} is not a known axis. Choices are: {sorted(all_axes.keys())}"
        )

nominal_cols = list(all_axes.keys()) if args.allaxes else list(args.axes)
if not nominal_cols:
    nominal_cols = [arm.mass_col]
logger.info("Histogramming %d axes.", len(nominal_cols))

# Columns that are per-candidate vectors need a scalar companion to fill a
# histogram. Everything the axis map keys on is either a candidate column or an
# already-scalar event quantity; the split is resolved per dataset in the graph.
cutflows = {}
cand_cutflows = {}
hist_names = set()


def _selection_list(df_is_data):
    """The ordered selection list.

    Sentinel rejection comes before every kinematic cut on the arm, because a
    sentinel passes a one-sided bound. The quality gate follows immediately,
    for the same reason in spirit: a fit that did not converge has no mass
    worth cutting on.
    """
    gate_bounds = {
        "edmRef": args.gateEdmRef,
        "chisq": args.gateChisq,
        "ndof": args.gateNdof,
        "massErr": args.gateMassErr,
    }
    dimuon_source = dimuon_vtx_arm(arm)

    selections = [
        (f"{arm.name} fit ok", lambda df: sel.select_arm_ok(df, arm)),
        (f"{arm.name} no sentinel", lambda df: sel.select_arm_finite(df, arm)),
    ]
    selections += sel.quality_gate_selections(arm, gate_bounds)
    selections += [
        ("opposite-sign muons", sel.select_opposite_sign_muons),
        (
            f"muon pT > {args.muonMinPt:g}",
            lambda df: sel.select_muon_pt(df, args.muonMinPt),
        ),
        (
            f"muon |eta| < {args.muonMaxEta:g}",
            lambda df: sel.select_muon_eta(df, args.muonMaxEta),
        ),
    ]
    if args.requireGlobalMuons:
        selections.append(("both muons global", sel.select_muon_global))
    if args.minDimuonPt > 0:
        selections.append(
            (
                f"dimuon pT > {args.minDimuonPt:g}",
                lambda df: sel.select_dimuon_pt(df, args.minDimuonPt),
            )
        )
    selections += [
        (
            f"kaon pT > {args.bachelorMinPt:g}",
            lambda df: sel.select_bachelor_pt(df, args.bachelorMinPt),
        ),
        (
            f"kaon |eta| < {args.bachelorMaxEta:g}",
            lambda df: sel.select_bachelor_eta(df, args.bachelorMaxEta),
        ),
        (
            f"mass in [{args.massLow:g}, {args.massHigh:g}]",
            lambda df: sel.select_mass_window(df, arm, args.massLow, args.massHigh),
        ),
        (
            f"nLegsRefit >= {args.minLegsRefit}",
            lambda df: sel.select_n_legs_refit(df, args.minLegsRefit),
        ),
    ]
    if args.minDimuonVtxProb > 0:
        selections.append(
            (
                f"dimuon vtx prob > {args.minDimuonVtxProb:g} ({dimuon_source.name})",
                lambda df: sel.select_dimuon_vertex_probability(
                    df, arm, args.minDimuonVtxProb
                ),
            )
        )
    if args.minSl3d > 0:
        selections.append(
            (
                f"{arm.name} sl3d > {args.minSl3d:g}",
                lambda df: sel.select_sl3d(df, arm, args.minSl3d),
            )
        )
    if args.maxAlphaBS >= 0:
        selections.append(
            (
                f"{arm.name} alphaBS < {args.maxAlphaBS:g}",
                lambda df: sel.select_alpha_bs(df, arm, args.maxAlphaBS),
            )
        )
    if args.minDimuonSl3d > 0:
        selections.append(
            (
                f"dimuon sl3d > {args.minDimuonSl3d:g} ({dimuon_source.name})",
                lambda df: sel.select_dimuon_sl3d(df, arm, args.minDimuonSl3d),
            )
        )
    if args.maxDimuonAlphaBS >= 0:
        selections.append(
            (
                f"dimuon alphaBS < {args.maxDimuonAlphaBS:g} ({dimuon_source.name})",
                lambda df: sel.select_dimuon_alpha_bs(df, arm, args.maxDimuonAlphaBS),
            )
        )
    if args.bachelorMaxPt > 0:
        selections.append(
            (
                f"kaon pT < {args.bachelorMaxPt:g}",
                lambda df: sel.select_bachelor_pt_max(df, args.bachelorMaxPt),
            )
        )
    # Candidate vertex probability. Declared since the first version of this
    # histmaker and never applied -- so every result so far was produced with it
    # silently off, including at the 0.3 the 2018 selection used, which is that
    # selection's tightest cut.
    if args.minVtxProb > 0:
        selections.append(
            (
                f"{arm.name} vtx prob > {args.minVtxProb:g}",
                lambda df: sel.select_vertex_probability(
                    df, f"cand_{arm.name}_vtxProb", args.minVtxProb
                ),
            )
        )
    # Bachelor and joint-fit quality, each its own cutflow line so it is visible
    # which one bites.
    if args.bachelorMinValidHits >= 0:
        selections.append(
            (
                f"kaon nValidHits > {args.bachelorMinValidHits}",
                lambda df: sel.select_bachelor_n_valid_hits(
                    df, args.bachelorMinValidHits
                ),
            )
        )
    if args.minBachelorAbsD0 >= 0:
        selections.append(
            (
                f"kaon |d0| > {args.minBachelorAbsD0:g}",
                lambda df: sel.select_bachelor_abs_d0(df, args.minBachelorAbsD0),
            )
        )
    if args.minBachelorDedx >= 0:
        selections.append(
            (
                f"kaon dE/dx > {args.minBachelorDedx:g}",
                lambda df: sel.select_bachelor_dedx(df, args.minBachelorDedx),
            )
        )
    if args.maxJointCvhNiter >= 0:
        selections.append(
            (
                f"{arm.name} nIter <= {args.maxJointCvhNiter}",
                lambda df: sel.select_joint_n_iter(df, arm, args.maxJointCvhNiter),
            )
        )
    if args.minBachelorAbsD0Sig > 0:
        selections.append(
            (
                f"bachelor |d0|/sigma > {args.minBachelorAbsD0Sig:g}",
                lambda df: sel.select_bachelor_abs_d0_sig(df, args.minBachelorAbsD0Sig),
            )
        )
    if args.maxChisqNdof > 0:
        selections.append(
            (
                f"{arm.name} chisq/ndof < {args.maxChisqNdof:g}",
                lambda df: sel.select_chisq_ndof(df, args.maxChisqNdof),
            )
        )
    return selections


# Pileup reweighting for the simulation, on the generated number of true
# interactions (gap G4, closed). What was previously reweighted -- the
# reconstructed vertex count -- is now left alone, so its data/MC agreement is
# a result rather than something imposed.
pileup_weights = None
if args.pileupWeights:
    with open(args.pileupWeights) as fh:
        pileup_weights = {int(k): float(v) for k, v in json.load(fh).items()}
    logger.info(
        "Loaded %d nTrueInt weights from %s (range %.3f - %.3f)",
        len(pileup_weights),
        args.pileupWeights,
        min(pileup_weights.values()),
        max(pileup_weights.values()),
    )


def build_graph(df, dataset):
    logger.info("build graph for dataset: %s", dataset.name)
    results = []
    cutflow = {}
    storage_type = hist.storage.Double()

    if dataset.is_data:
        df = df.DefinePerSample("weight", "1.0")
    else:
        # GenEvt_weight, not genWeight: this NanoAOD renames it. It is
        # identically 1.0 in this production, so the sum of weights is an
        # event count -- see the normalisation section of the proposal.
        df = df.Define("weight", "GenEvt_weight")
    df = sel.define_pileup_weight(df, dataset.is_data, pileup_weights)
    df = df.Define("nominal_weight", "static_cast<double>(weight) * cand_pileupWeight")

    weightsum = df.SumAndCount("weight")
    cutflow["Total"] = weightsum[0]

    # The gen-level denominator, booked HERE and not later. Every filter below
    # this line -- trigger, candidate multiplicity, kinematics -- is part of
    # what the efficiency is measuring, so a denominator defined after any of
    # them measures nothing.
    if args.genDiagnostics and not dataset.is_data:
        df = sel.define_gen_truth_signal(df, dataset.is_data)
        results.append(
            df.HistoBoost(
                "genTruthKaon",
                [
                    axis_gen_truth_pt,
                    axis_gen_truth_eta,
                    hist.axis.Boolean(name="inAcc"),
                ],
                [
                    "genTruth_kaonPt",
                    "genTruth_kaonEta",
                    "genTruth_inAcc",
                    "nominal_weight",
                ],
                storage=hist.storage.Double(),
            )
        )
        hist_names.add("genTruthKaon")

    df = sel.apply_trigger(df, trigger_column)
    cutflow[trigger_column] = df.SumAndCount("weight")[0]

    df = df.Filter(f"n{CAND} > 0", "at least one candidate")
    cutflow["candidate > 0"] = df.SumAndCount("weight")[0]

    df = sel.init_candidate_mask(df)
    df = sel.define_leg_kinematics(df)
    df = sel.define_raw_dimuon(df)

    # Both of these are defined before the selections, not after, because they are
    # now cut on as well as ranked on. The vertex probability used to be computed
    # only for the candidate ranking, which is why --minVtxProb existed for months
    # without doing anything.
    df, prob_col = sel.define_arm_vertex_probability(df, arm)
    df = sel.define_chisq_ndof(df, arm)

    cand_cutflow = {}
    df, cand_start = sel.count_candidates(df, "start")
    cand_cutflow["candidate > 0"] = cand_start

    df, cutflow, dfs_per_cut = sel.apply_selections(
        df, _selection_list(dataset.is_data), cutflow, cand_cutflow
    )

    # Rank candidates by the nominal arm's vertex probability.
    df = sel.select_best_candidate(df, prob_col)
    cutflow["best candidate"] = df.SumAndCount("weight")[0]

    df = sel.define_gen_category(df, dataset.is_data)
    df = sel.define_bachelor_pt_response(df, dataset.is_data)

    # Undo the generator's forced b-hadron decays, so the stacked components
    # carry their physical relative rates. This changes composition, never a
    # single component's shape, and it is applied after the best-candidate
    # choice because it is keyed to that candidate's matched species.
    df = sel.define_forced_decay_weight(
        df, dataset.is_data, not args.noForcedDecayWeight
    )
    df = df.Redefine(
        "nominal_weight",
        "static_cast<double>(weight) * cand_pileupWeight * cand_forcedDecayWeight",
    )

    # Scalar companions for every candidate-level column being histogrammed.
    fill_cols = {}
    for var in nominal_cols:
        if var in MC_ONLY_AXES and dataset.is_data:
            continue
        if var in EVENT_SCALAR_AXES:
            fill_cols[var] = var
            continue
        if var == "PV_z":
            # PV_z is a per-vertex vector; the leading vertex is the one the
            # candidate geometry is referenced to.
            df = df.Define("PV_z_scalar", "PV_z.size() ? PV_z[0] : -999.")
            fill_cols[var] = "PV_z_scalar"
            continue
        if not df.HasColumn(var):
            logger.warning(
                "Column %s not available for dataset %s, skipping", var, dataset.name
            )
            continue
        coltype = df.GetColumnType(var)
        if "RVec" in coltype or "vector<" in coltype:
            alias = f"{var}_scalar"
            df = sel.define_scalar(df, var, alias)
            fill_cols[var] = alias
        else:
            fill_cols[var] = var

    for var, col in fill_cols.items():
        name = f"nominal_{var}"
        results.append(
            df.HistoBoost(
                name,
                [all_axes[var], axis_gen_category],
                [col, "cand_genCategory", "nominal_weight"],
                storage=storage_type,
            )
        )
        hist_names.add(name)

    # -----------------------------------------------------------------------
    # A/e/M momentum-scale variations
    # -----------------------------------------------------------------------
    if scale_unc_helper is not None and not dataset.is_data:
        df, scale_prefix = sel.define_scale_variation_inputs(
            df,
            dataset.is_data,
            legs=(
                sel.SCALE_VAR_LEGS
                if args.scaleVarLegs == "all"
                else sel.SCALE_VAR_LEGS[:1]
            ),
            muon_source=args.scaleVarMuonSource,
        )
        df, scale_cols = muon_calibration.muon_reweight_helper_cols(
            df, scale_unc_helper, scale_prefix, None
        )
        df = df.Define(
            "nominal_muonScaleSyst_responseWeights_tensor",
            scale_unc_helper,
            [*scale_cols, "nominal_weight"],
        )
        # Booked over the fit binning, not over mass alone: the deliverable is
        # the shift the variation induces per fit cell, which is what says
        # whether the channel has constraining power.
        for col in (f"{CAND}_kaonPt", f"{CAND}_kaonEta", f"{CAND}_charge"):
            if not df.HasColumn(f"{col}_scalar"):
                df = sel.define_scalar(df, col, f"{col}_scalar")
        results.append(
            df.HistoBoost(
                "nominal_muonScaleSyst_responseWeights",
                [
                    axis_fit_bachelor_pt,
                    axis_fit_bachelor_eta,
                    axis_fit_bachelor_charge,
                    all_axes[arm.mass_col],
                    axis_gen_category,
                ],
                [
                    f"{CAND}_kaonPt_scalar",
                    f"{CAND}_kaonEta_scalar",
                    f"{CAND}_charge_scalar",
                    f"{arm.mass_col}_scalar",
                    "cand_genCategory",
                    "nominal_muonScaleSyst_responseWeights_tensor",
                ],
                tensor_axes=scale_unc_helper.tensor_axes,
                # Weight, not Double: the tensor writer pushes this histogram
                # through the *same* reduction path as `fitbins`
                # (rabbit_cvhnano_helpers.reduce_fitbins), and that path folds
                # the pT overflow and rebins by operating on values and
                # variances together. Matching storage is what lets the
                # variation take a byte-identical route to the nominal instead
                # of a parallel implementation that can drift from it.
                storage=hist.storage.Weight(),
            )
        )
        hist_names.add("nominal_muonScaleSyst_responseWeights")

        # The gen-match audit, counted on whatever sample is actually run rather
        # than asserted from 12 files. Every candidate in the exclusive-signal
        # and other-real-b categories was measured to have all three legs
        # matched, and no not-fully-matched candidate was signal; this is the
        # histogram that says so on the full statistics.
        #
        # A hard throw would have been the natural assertion, but a throw inside
        # an RDataFrame Define kills the event loop -- the same reason the exact
        # shift is guarded. So the condition is counted and checked after the
        # run instead.
        results.append(
            df.HistoBoost(
                "scaleVarMatchAudit",
                [axis_gen_category, hist.axis.Boolean(name="allMatched")],
                [
                    "cand_genCategory",
                    f"{scale_prefix}_allMatched",
                    "nominal_weight",
                ],
                storage=hist.storage.Double(),
            )
        )
        hist_names.add("scaleVarMatchAudit")

    # -----------------------------------------------------------------------
    # Generator-truth diagnostics on the reconstructed candidates
    # -----------------------------------------------------------------------
    if args.genDiagnostics and not dataset.is_data:
        df = sel.define_gen_match_diagnostics(df, dataset.is_data)
        # The scale-variation block defines this too, but the diagnostics must
        # not depend on --scaleVariations being on to work.
        if not df.HasColumn(f"{CAND}_kaonEta_scalar"):
            df = sel.define_scalar(df, f"{CAND}_kaonEta", f"{CAND}_kaonEta_scalar")
        # The mass distribution split by how many legs the matcher accepted. The
        # no-ancestor ("combinatorial") category is filled both by genuinely
        # random track triplets and by real decays one of whose legs the matcher
        # rejected, and those two are not the same background: the second peaks
        # at the B mass. Only this histogram separates them.
        results.append(
            df.HistoBoost(
                "genMatchDiag",
                [
                    all_axes[arm.mass_col],
                    axis_gen_category,
                    axis_n_legs_matched,
                    axis_ancestor_state,
                ],
                [
                    f"{arm.mass_col}_scalar",
                    "cand_genCategory",
                    "cand_nLegsMatched",
                    "cand_ancestorState",
                    "nominal_weight",
                ],
                storage=hist.storage.Double(),
            )
        )
        hist_names.add("genMatchDiag")

        # The quantity the matcher cut on. A hard edge at +-0.1 with population
        # piled against it means the cut is removing real candidates rather than
        # separating right from wrong.
        results.append(
            df.HistoBoost(
                "bachelorRelDPt",
                [
                    axis_rel_dpt,
                    axis_gen_truth_pt,
                    axis_fit_bachelor_eta,
                    axis_gen_category,
                ],
                [
                    "cand_bachelorRelDPt",
                    "cand_bachelorGenPt",
                    f"{CAND}_kaonEta_scalar",
                    "cand_genCategory",
                    "nominal_weight",
                ],
                storage=hist.storage.Double(),
            )
        )
        hist_names.add("bachelorRelDPt")

        # The reconstructed bachelor spectrum of candidates classified signal,
        # against the gen denominator above, in the same binning.
        results.append(
            df.HistoBoost(
                "genTruthKaonReco",
                [axis_gen_truth_pt, axis_gen_truth_eta, axis_gen_category],
                [
                    "cand_bachelorGenPt",
                    f"{CAND}_kaonEta_scalar",
                    "cand_genCategory",
                    "nominal_weight",
                ],
                storage=hist.storage.Double(),
            )
        )
        hist_names.add("genTruthKaonReco")

    # A bare "nominal" so plotting tools that expect one find it.
    results.append(
        df.HistoBoost(
            "nominal",
            [all_axes[arm.mass_col], axis_gen_category],
            [f"{arm.mass_col}_scalar", "cand_genCategory", "nominal_weight"],
            storage=storage_type,
        )
    )

    # The differential fit binning: mass in bins of bachelor (pT, eta, charge).
    # Booked on fine, regular pT and mass axes rather than the eventual quantile
    # binning, because quantile edges have to be derived *from* this histogram.
    # eta is on the intended 0.1-wide binning already, since that one is fixed.
    if args.qualityHists:
        for col, axis in ID_AXES.items():
            if not df.HasColumn(f"{col}_scalar"):
                df = sel.define_scalar(df, col, f"{col}_scalar")
            if not df.HasColumn(f"{CAND}_kaonPt_scalar"):
                df = sel.define_scalar(df, f"{CAND}_kaonPt", f"{CAND}_kaonPt_scalar")
            results.append(
                df.HistoBoost(
                    f"quality_{axis.name}",
                    [axis, axis_id_bachelor_pt, axis_gen_category],
                    [
                        f"{col}_scalar",
                        f"{CAND}_kaonPt_scalar",
                        "cand_genCategory",
                        "nominal_weight",
                    ],
                    storage=hist.storage.Weight(),
                )
            )

    if args.fitBinsHist:
        # Defined here rather than relying on the --allaxes loop, so the
        # histogram does not silently depend on which axes were requested.
        for col in (f"{CAND}_kaonPt", f"{CAND}_kaonEta", f"{CAND}_charge"):
            if not df.HasColumn(f"{col}_scalar"):
                df = sel.define_scalar(df, col, f"{col}_scalar")
        results.append(
            df.HistoBoost(
                "fitbins",
                [
                    axis_fit_bachelor_pt,
                    axis_fit_bachelor_eta,
                    axis_fit_bachelor_charge,
                    all_axes[arm.mass_col],
                    axis_gen_category,
                ],
                [
                    f"{CAND}_kaonPt_scalar",
                    f"{CAND}_kaonEta_scalar",
                    f"{CAND}_charge_scalar",
                    f"{arm.mass_col}_scalar",
                    "cand_genCategory",
                    "nominal_weight",
                ],
                # Weight, not Double: the tensor stage needs variances to
                # separate the simulated template's statistical contribution
                # from the data's, and the simulation is genuinely weighted --
                # the forced-decay weight spans 0.048 to 51.8 across species,
                # so variances are not counts. Only this histogram changes;
                # everything else stays on the shared Double storage.
                storage=hist.storage.Weight(),
            )
        )

    cutflows[dataset.name] = cutflow
    cand_cutflows[dataset.name] = cand_cutflow
    return results, weightsum


resultdict = narf.build_and_run(datasets, build_graph)

for dataset_name, actions in cutflows.items():
    resultdict[dataset_name]["cutflow"] = {
        name: action.GetValue() for name, action in actions.items()
    }
for dataset_name, actions in cand_cutflows.items():
    resultdict[dataset_name]["cand_cutflow"] = {
        name: action.GetValue() for name, action in actions.items()
    }

# -------------------------------------------------------------------------
# Post-run audits of the scale variations
# -------------------------------------------------------------------------
if scale_unc_helper is not None:
    n_fallback = scale_unc_helper.nExactFallback()
    if n_fallback:
        logger.warning(
            "The exact mass-aware shift was unphysical %d times and fell back "
            "to the massless linearisation for that (leg, variation). Above "
            "1 GeV this is expected to be zero; a non-zero count means the "
            "bachelor pT floor is low enough for E - m_K to be comparable to "
            "the e uncertainty.",
            n_fallback,
        )
    else:
        logger.info("Exact mass-aware shift: no unphysical evaluations.")

    for dataset_name, payload in resultdict.items():
        audit = payload.get("output", {}).get("scaleVarMatchAudit")
        if audit is None:
            continue
        h = audit.get() if hasattr(audit, "get") else audit
        values = h.values()
        # allMatched is the last axis: index 0 is False, 1 is True.
        for icat in range(values.shape[0]):
            unmatched, matched = values[icat, 0], values[icat, 1]
            total = unmatched + matched
            if total <= 0:
                continue
            logger.info(
                "gen-match audit, %s, genCategory %d: %.0f of %.0f "
                "candidates fully matched (%.2f%%)",
                dataset_name,
                icat,
                matched,
                total,
                100.0 * matched / total,
            )
        # The claim under test: the signal category cannot contain a partially
        # matched candidate, because the chain walk needs three matched legs to
        # reach a common ancestor.
        if values[GEN_SIGNAL, 0] > 0:
            logger.error(
                "%s: %.0f exclusive-signal candidates are NOT fully gen "
                "matched. That should be structurally impossible, so the truth "
                "categorisation and the leg matching disagree and the variation "
                "weights on those candidates are not meaningful.",
                dataset_name,
                values[GEN_SIGNAL, 0],
            )

if not args.noScaleToData:
    scale_to_data(resultdict)
    aggregate_groups(datasets, resultdict, args.aggregateGroups)

write_analysis_output(
    resultdict, f"{os.path.basename(__file__).replace('py', 'hdf5')}", args
)


def render_cutflow(key, title):
    """Render one cutflow table as text, or None if nothing was recorded."""
    samples = [n for n in resultdict if key in resultdict[n]]
    if not samples:
        return None
    names = []
    for sample in samples:
        for name in resultdict[sample][key]:
            if name not in names:
                names.append(name)

    header = f"{'Selection':<42}" + "".join(f"{s[:22]:>24}" for s in samples)
    lines = [title, "", header, "-" * len(header)]
    first = {s: None for s in samples}
    for name in names:
        row = f"{name:<42}"
        for sample in samples:
            val = resultdict[sample][key].get(name)
            if val is None:
                row += f"{'-':>24}"
                continue
            if first[sample] is None:
                first[sample] = val
            frac = val / first[sample] if first[sample] else 0.0
            row += f"{val:>15.4g} ({frac:5.3f})"
        lines.append(row)
    return "\n".join(lines)


if args.cutflow or args.saveCutflow:
    tables = [
        render_cutflow(
            "cutflow",
            "Cutflow, EVENTS with at least one surviving candidate "
            "(yield, fraction of total):",
        ),
        render_cutflow(
            "cand_cutflow",
            "Cutflow, CANDIDATES surviving (yield, fraction of candidates "
            "entering). An event with several candidates survives a cut that "
            "kills most of them, so per-candidate cuts bite harder here than "
            "in the event table above.",
        ),
    ]
    table = "\n\n".join(t for t in tables if t)
    print()
    print(table)
    print()

    if args.saveCutflow:
        os.makedirs(args.saveCutflow, exist_ok=True)
        postfix = args.postfix or "cutflow"
        path = os.path.join(args.saveCutflow, f"cutflow_{postfix}.txt")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(table + "\n")
        logger.info("Cutflow written to %s", path)

print(
    "hist variable names for plotting:\n\n "
    + " ".join(sorted(n.replace("nominal_", "") for n in hist_names))
    + "\n"
)
