import os

import hist
import ROOT

import narf
from wremnants.production import muon_calibration
from wremnants.production.histmaker_tools import write_analysis_output
from wremnants.utilities import common, parsing
from wums import logging

analysis_label = common.analysis_label(os.path.basename(__file__))
parser, initargs = parsing.common_parser(analysis_label)
parser.add_argument(
    "--resonance",
    type=str,
    choices=["jpsi", "upsilon"],
    help="Resonance for selection",
    default="upsilon",
)
parser.add_argument(
    "--etaBins",
    type=int,
    default=None,
    help="Override the number of bins for each reconstructed muon eta axis",
)
parser.add_argument(
    "--ptBins",
    type=int,
    choices=[2, 3],
    default=None,
    help=(
        "Override the number of bins for each reconstructed muon pt axis. "
        "Currently 2 and 3 are implemented for J/psi fit studies."
    ),
)
parser.add_argument(
    "--pixelHitSelection",
    choices=["all", "bothPixelHits", "anyMissingPixelHits"],
    default="all",
    help=(
        "Optional selection on the number of valid pixel hits in the track refit. "
        "Use bothPixelHits for both muons passing --minValidPixelHits; "
        "use anyMissingPixelHits for the complementary category."
    ),
)
parser.add_argument(
    "--minValidPixelHits",
    type=int,
    default=1,
    help=(
        "Minimum valid pixel hits required per muon when "
        "--pixelHitSelection bothPixelHits is used."
    ),
)
parser.add_argument(
    "--deltaPhiSign",
    choices=["all", "positive", "negative"],
    default="all",
    help=(
        "Optional selection on the signed dimuon azimuthal separation, defined "
        "as TVector2::Phi_mpi_pi(phi_plus - phi_minus). Use negative for "
        "DeltaPhi < 0 and positive for DeltaPhi > 0."
    ),
)
parser.add_argument(
    "--resolutionPrefitUncertainty",
    type=float,
    default=3.0,
    help=(
        "Relative prefit uncertainty assigned to each fitted muon-resolution "
        "parameter when --fitMuonScaleAndResolution is used."
    ),
)
parser.add_argument(
    "--resolutionPrefitUncertaintyA",
    type=float,
    default=None,
    help=(
        "Override the relative prefit uncertainty for only the fitted "
        "muon-resolution a parameter. If omitted, --resolutionPrefitUncertainty "
        "is used for a, b, c, and d."
    ),
)
parser = parsing.set_parser_default(parser, "theoryCorr", [])
parser = parsing.set_parser_default(parser, "scale_A", 5.0)
parser = parsing.set_parser_default(parser, "scale_e", 5.0)

args = parser.parse_args()

if args.etaBins is not None and args.etaBins <= 0:
    raise ValueError("--etaBins must be a positive integer")
if args.ptBins is not None and args.ptBins <= 0:
    raise ValueError("--ptBins must be a positive integer")
if args.minValidPixelHits < 1:
    raise ValueError("--minValidPixelHits must be at least 1")
if args.resolutionPrefitUncertainty <= 0:
    raise ValueError("--resolutionPrefitUncertainty must be positive")
if (
    args.resolutionPrefitUncertaintyA is not None
    and args.resolutionPrefitUncertaintyA <= 0
):
    raise ValueError("--resolutionPrefitUncertaintyA must be positive")
if args.etaBins is not None and not args.fitMuonScaleAndResolution:
    raise ValueError(
        "--etaBins currently requires --fitMuonScaleAndResolution so scale "
        "uncertainties can be correlated across the requested eta groups"
    )

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

if args.muonScaleVariation not in ["onnxReweight", "smearingWeightsSplines"]:
    raise ValueError(
        "dimuon_resonances_calinput.py currently supports "
        "--muonScaleVariation onnxReweight or smearingWeightsSplines"
    )

calib_filepaths = common.calib_filepaths
scale_diff_weights_helper = (
    ROOT.wrem.SplinesDifferentialWeightsHelper(calib_filepaths["tflite_file"])
    if args.muonScaleVariation == "smearingWeightsSplines"
    else None
)
resolution_diff_weights_helper = (
    ROOT.wrem.SplinesDifferentialWeightsHelper(
        calib_filepaths["tflite_file_nosmearing"]
    )
    if args.muonScaleVariation == "smearingWeightsSplines"
    else None
)
(
    _mc_jpsi_crctn_helper,
    _data_jpsi_crctn_helper,
    _mc_jpsi_crctn_unc_helper,
    data_jpsi_crctn_unc_helper,
) = muon_calibration.make_jpsi_crctn_helpers(
    calib_filepaths,
    muon_corr_mc=args.muonCorrMC,
    muon_corr_data=args.muonCorrData,
    scale_var_method=args.muonScaleVariation,
    scale_A=args.scale_A,
    scale_e=args.scale_e,
    scale_M=args.scale_M,
    make_uncertainty_helper=True,
    include_covariance=not args.fitMuonScaleAndResolution,
    smearing=not args.noSmearing,
    fit_muon_scale=args.fitMuonScaleAndResolution,
    variation_eta_bins=args.etaBins,
)

if data_jpsi_crctn_unc_helper is None:
    raise ValueError(
        "Muon scale systematics require --muonCorrData massfit or lbl_massfit"
    )

_smearing_helper, smearing_uncertainty_helper = (
    (None, None)
    if args.noSmearing
    else muon_calibration.make_muon_smearing_helpers(
        scale_var_method=args.muonScaleVariation,
        parameter_variations=True,
        fit_muon_resolution=args.fitMuonScaleAndResolution,
        variation_eta_bins=args.etaBins,
        resolution_prefit_uncertainties=[
            args.resolutionPrefitUncertaintyA
            or args.resolutionPrefitUncertainty,
            args.resolutionPrefitUncertainty,
            args.resolutionPrefitUncertainty,
            args.resolutionPrefitUncertainty,
        ],
        resolution_prefit_uncertainties_mode="relative",
    )
)

(
    _pixel_multiplicity_helper,
    pixel_multiplicity_uncertainty_helper,
    pixel_multiplicity_uncertainty_helper_stat,
) = muon_calibration.make_pixel_multiplicity_helpers()


local_resonance_files = {
    "jpsi": {
        "data": ["/scratch/submit/cms/emanca/jpsicor_data.root"],
        "mc": [
            "/scratch/submit/cms/emanca/jpsicor_mc.root",
            "/scratch/submit/cms/emanca/jcor_mc_0to8.root",
        ],
    },
    "upsilon": {
        "data": ["/scratch/submit/cms/emanca/upsilon_data.root"],
        "mc": [
            "/scratch/submit/cms/emanca/upsilon_0to8.root",
            "/scratch/submit/cms/emanca/upsilon_8to13.root",
            "/scratch/submit/cms/emanca/upsilon_13toInf.root",
        ],
    },
}

trigger_channels = {
    "jpsi": [
        {
            "label": "dimuon20_jpsi",
            "cut": "HLT_Dimuon20_Jpsi",
            "mc": local_resonance_files["jpsi"]["mc"],
            "layer_corrected": True,
        },
        {
            "label": "doublemu4_jpsitrk_displaced",
            "cut": "HLT_DoubleMu4_JpsiTrk_Displaced",
            "mc": ["/scratch/submit/cms/emanca/BuToJpsiK_BMuonFilter_v2_BPH.root"],
            "layer_corrected": False,
        },
    ],
    "upsilon": [
        {
            "label": "inclusive",
            "cut": "1.",
            "mc": local_resonance_files["upsilon"]["mc"],
            "layer_corrected": True,
        },
    ],
}


def limited_files(files):
    if args.maxFiles is not None and args.maxFiles > 0:
        return files[: args.maxFiles]
    return files


def bool_filter(expression):
    return f"static_cast<bool>({expression})"


datasets = []
dataset_channels = {}
for channel in trigger_channels[args.resonance]:
    for sample_type in ["data", "mc"]:
        dataset_name = f"{args.resonance}_{sample_type}_{channel['label']}"
        dataset_channels[dataset_name] = channel
        datasets.append(
            narf.Dataset(
                name=dataset_name,
                filepaths=limited_files(
                    local_resonance_files[args.resonance]["data"]
                    if sample_type == "data"
                    else channel["mc"]
                ),
                is_data=sample_type == "data",
                xsec=None if sample_type == "data" else 1.0,
                group="Data" if sample_type == "data" else args.resonance.upper(),
            )
        )


resonance_options = {
    "jpsi": {
        "name": "JPsi",
        "default_eta_bins": 24,
        "eta_range": (-2.4, 2.4),
        "pt_axis": hist.axis.Variable([4.2, 7.0, 10.5, 15.0, 25.0], name="pt1"),
        "pt2_axis": hist.axis.Variable([4.2, 7.0, 10.5, 15.0, 25.0], name="pt2"),
        "mass_axis": hist.axis.Variable(
            [
                2.92,
                2.965,
                3.000,
                3.030,
                3.050,
                3.065,
                3.077,
                3.087,
                3.091,
                3.095,
                3.101,
                3.106,
                3.111,
                3.116,
                3.121,
                3.126,
                3.132,
                3.139,
                3.148,
                3.159,
                3.173,
                3.190,
                3.210,
                3.233,
                3.258,
                3.28,
            ],
            name="mass",
        ),
        "mass_range": (2.8, 3.35),
    },
    "upsilon": {
        "name": "Y",
        "default_eta_bins": 8,
        "eta_range": (-0.8, 0.8),
        "pt_axis": hist.axis.Variable([4.2, 6.0, 7.9, 10.3, 25.0], name="pt1"),
        "pt2_axis": hist.axis.Variable([4.2, 6.0, 7.9, 10.3, 25.0], name="pt2"),
        "mass_axis": hist.axis.Regular(25, 9.0, 9.7, name="mass"),
        "mass_range": (8.8, 9.6),
    },
}

cfg = resonance_options[args.resonance]
eta_bins = args.etaBins or cfg["default_eta_bins"]
eta_min, eta_max = cfg["eta_range"]
mass_axis = cfg["mass_axis"]

if args.ptBins == 2:
    if args.resonance == "jpsi":
        pt_axis = hist.axis.Variable([4.2, 10.5, 25.0], name="pt1")
        pt2_axis = hist.axis.Variable([4.2, 10.5, 25.0], name="pt2")
    else:
        pt_axis = hist.axis.Variable([4.2, 7.9, 25.0], name="pt1")
        pt2_axis = hist.axis.Variable([4.2, 7.9, 25.0], name="pt2")
elif args.ptBins == 3:
    if args.resonance == "jpsi":
        pt_axis = hist.axis.Variable([4.2, 7.0, 10.5, 25.0], name="pt1")
        pt2_axis = hist.axis.Variable([4.2, 7.0, 10.5, 25.0], name="pt2")
    else:
        pt_axis = hist.axis.Variable([4.2, 7.9, 10.3, 25.0], name="pt1")
        pt2_axis = hist.axis.Variable([4.2, 7.9, 10.3, 25.0], name="pt2")
else:
    pt_axis = cfg["pt_axis"]
    pt2_axis = cfg["pt2_axis"]

calibration_axes = [
    hist.axis.Regular(eta_bins, eta_min, eta_max, name="eta1"),
    hist.axis.Regular(eta_bins, eta_min, eta_max, name="eta2"),
    pt_axis,
    pt2_axis,
    mass_axis,
]


def reco_columns(channel):
    suffix = "cor" if channel["layer_corrected"] else ""
    return {
        "plus_pt": f"Muplus{suffix}_pt",
        "minus_pt": f"Muminus{suffix}_pt",
        "plus_eta": f"Muplus{suffix}_eta",
        "minus_eta": f"Muminus{suffix}_eta",
        "plus_phi": f"Muplus{suffix}_phi",
        "minus_phi": f"Muminus{suffix}_phi",
        "mass": f"Jpsi{suffix}_mass",
    }


def reco_selection(cols):
    mass_min, mass_max = cfg["mass_range"]
    return (
        f"{cols['plus_pt']} > 1.0 && {cols['minus_pt']} > 1.0 && "
        f"{cols['plus_pt']} < 100.0 && {cols['minus_pt']} < 100.0 && "
        f"std::fabs({cols['plus_eta']}) < {eta_max} && "
        f"std::fabs({cols['minus_eta']}) < {eta_max} && "
        f"{cols['mass']} > {mass_min} && {cols['mass']} < {mass_max}"
    )


def pixel_hit_selection():
    if args.pixelHitSelection == "all":
        return "1."
    if args.pixelHitSelection == "bothPixelHits":
        return (
            f"Muplus_nvalidpixel >= {args.minValidPixelHits} && "
            f"Muminus_nvalidpixel >= {args.minValidPixelHits}"
        )
    if args.pixelHitSelection == "anyMissingPixelHits":
        return (
            f"Muplus_nvalidpixel < {args.minValidPixelHits} || "
            f"Muminus_nvalidpixel < {args.minValidPixelHits}"
        )
    raise ValueError(f"Unsupported --pixelHitSelection {args.pixelHitSelection}")


def delta_phi_selection(cols):
    if args.deltaPhiSign == "all":
        return "1."

    comparison = "<" if args.deltaPhiSign == "negative" else ">"
    return (
        f"TVector2::Phi_mpi_pi({cols['plus_phi']} - {cols['minus_phi']}) "
        f"{comparison} 0.0"
    )


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")

    results = []
    channel = dataset_channels[dataset.name]
    reco_cols = reco_columns(channel)
    calibration_cols = [
        reco_cols["plus_eta"],
        reco_cols["minus_eta"],
        reco_cols["plus_pt"],
        reco_cols["minus_pt"],
        reco_cols["mass"],
    ]

    df = df.DefinePerSample("weight", "1.0")
    weightsum = df.SumAndCount("weight")

    df = df.Filter(reco_selection(reco_cols))
    df = df.Filter(bool_filter(pixel_hit_selection()))
    df = df.Filter(bool_filter(delta_phi_selection(reco_cols)))
    if not dataset.is_data:
        df = df.Filter(
            "Jpsigen_mass > 0.0 && Muplusgen_pt > 0.0 && Muminusgen_pt > 0.0"
        )
    df = df.Filter(bool_filter(channel["cut"]))

    pixel_multiplicity_cols = [
        "pixel_triggerCat",
        "pixel_eta",
        "pixel_pt",
        "pixel_charge",
        "pixel_nvalidpixel",
    ]
    if not dataset.is_data:
        df = (
            df
            # The calInput ntuples have no per-leg trigger matching; both selected
            # paths are dimuon triggers, so both candidates use the triggering map.
            .DefinePerSample(
                "pixel_triggerCat",
                "ROOT::VecOps::RVec<wrem::TriggerCat>{wrem::TriggerCat::triggering, wrem::TriggerCat::triggering}",
            )
            .Define(
                "pixel_eta",
                f"ROOT::VecOps::RVec<float>{{float({reco_cols['plus_eta']}), float({reco_cols['minus_eta']})}}",
            )
            .Define(
                "pixel_pt",
                f"ROOT::VecOps::RVec<float>{{float({reco_cols['plus_pt']}), float({reco_cols['minus_pt']})}}",
            )
            .Define("pixel_charge", "ROOT::VecOps::RVec<int>{1, -1}")
            .Define(
                "pixel_nvalidpixel",
                "ROOT::VecOps::RVec<int>{int(Muplus_nvalidpixel), int(Muminus_nvalidpixel)}",
            )
        )

    hist_name = f"{cfg['name']}_{'data' if dataset.is_data else 'mc'}"
    results.append(
        df.HistoBoost(
            hist_name,
            calibration_axes,
            [*calibration_cols, "weight"],
        )
    )

    if not dataset.is_data:
        df = (
            df.Define("nominal_weight", "weight")
            .Define(
                "scale_recoPt",
                f"ROOT::VecOps::RVec<float>{{float({reco_cols['plus_pt']}), float({reco_cols['minus_pt']})}}",
            )
            .Define(
                "scale_recoEta",
                f"ROOT::VecOps::RVec<float>{{float({reco_cols['plus_eta']}), float({reco_cols['minus_eta']})}}",
            )
            .Define(
                "scale_recoPhi",
                f"ROOT::VecOps::RVec<float>{{float({reco_cols['plus_phi']}), float({reco_cols['minus_phi']})}}",
            )
            .Define("scale_recoCharge", "ROOT::VecOps::RVec<int>{1, -1}")
            .Define(
                "scale_genPt",
                "ROOT::VecOps::RVec<float>{float(Muplusgen_pt), float(Muminusgen_pt)}",
            )
            .Define(
                "scale_genEta",
                "ROOT::VecOps::RVec<float>{float(Muplusgen_eta), float(Muminusgen_eta)}",
            )
            .Define(
                "scale_genPhi",
                "ROOT::VecOps::RVec<float>{float(Muplusgen_phi), float(Muminusgen_phi)}",
            )
            .Define("scale_genCharge", "ROOT::VecOps::RVec<int>{1, -1}")
            .Define("scale_muon_source", "ROOT::VecOps::RVec<int>{443, 443}")
        )

        input_kinematics = [
            "scale_recoPt",
            "scale_recoEta",
            "scale_recoCharge",
            "scale_genPt",
            "scale_genEta",
            "scale_genCharge",
        ]
        if scale_diff_weights_helper is not None:
            df = df.Define(
                "scale_response_weight",
                scale_diff_weights_helper,
                input_kinematics,
            )
        if resolution_diff_weights_helper is not None:
            df = df.Define(
                "scale_resolution_response_weight",
                resolution_diff_weights_helper,
                input_kinematics,
            )

        df, scale_cols = muon_calibration.muon_reweight_helper_cols(
            df,
            data_jpsi_crctn_unc_helper,
            "scale",
            "scale_response_weight",
        )
        df = df.Define(
            "nominal_muonScaleSyst_responseWeights_tensor",
            data_jpsi_crctn_unc_helper,
            [*scale_cols, "nominal_weight"],
        )
        results.append(
            df.HistoBoost(
                "nominal_muonScaleSyst_responseWeights",
                calibration_axes,
                [*calibration_cols, "nominal_muonScaleSyst_responseWeights_tensor"],
                tensor_axes=data_jpsi_crctn_unc_helper.tensor_axes,
                storage=hist.storage.Double(),
            )
        )

        df = muon_calibration.add_resolution_uncertainty(
            df,
            calibration_axes,
            results,
            calibration_cols,
            smearing_uncertainty_helper,
            "scale",
            storage_type=hist.storage.Double(),
            response_weight_col="scale_resolution_response_weight",
        )

        df = df.Define(
            "nominal_pixelMultiplicitySyst_tensor",
            pixel_multiplicity_uncertainty_helper,
            [*pixel_multiplicity_cols, "nominal_weight"],
        )
        results.append(
            df.HistoBoost(
                "nominal_pixelMultiplicitySyst",
                calibration_axes,
                [*calibration_cols, "nominal_pixelMultiplicitySyst_tensor"],
                tensor_axes=pixel_multiplicity_uncertainty_helper.tensor_axes,
                storage=hist.storage.Double(),
            )
        )

        if args.pixelMultiplicityStat:
            df = df.Define(
                "nominal_pixelMultiplicityStat_tensor",
                pixel_multiplicity_uncertainty_helper_stat,
                [*pixel_multiplicity_cols, "nominal_weight"],
            )
            results.append(
                df.HistoBoost(
                    "nominal_pixelMultiplicityStat",
                    calibration_axes,
                    [*calibration_cols, "nominal_pixelMultiplicityStat_tensor"],
                    tensor_axes=pixel_multiplicity_uncertainty_helper_stat.tensor_axes,
                    storage=hist.storage.Double(),
                )
            )

    return results, weightsum


resultdict = narf.build_and_run(datasets, build_graph, event_tree="tree")

fout = f"{os.path.basename(__file__).replace('py', 'hdf5')}"
name_append = [args.resonance, args.era]
if args.etaBins is not None:
    name_append.append(f"etaBins_{args.etaBins}")
if args.ptBins is not None:
    name_append.append(f"ptBins_{args.ptBins}")
if args.deltaPhiSign != "all":
    name_append.append(f"deltaPhi_{args.deltaPhiSign}")
write_analysis_output(resultdict, fout, args, name_append=name_append)
