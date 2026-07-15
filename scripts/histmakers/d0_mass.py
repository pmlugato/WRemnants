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
    "--dataFile",
    default="/scratch/submit/cms/emanca/D0data.root",
    help="Input ROOT file for data",
)
parser.add_argument(
    "--mcFile",
    default="/scratch/submit/cms/emanca/D0mc.root",
    help="Input ROOT file for MC",
)
parser.add_argument(
    "--treeName",
    default="tree",
    help="Input tree name",
)
parser.add_argument(
    "--etaBins",
    type=int,
    default=None,
    help="Override the number of eta bins used by fitted scale/resolution variations",
)
parser.add_argument(
    "--resolutionPrefitUncertainty",
    type=float,
    default=3.0,
    help=(
        "Relative prefit uncertainty assigned to each fitted resolution "
        "parameter when --fitMuonScaleAndResolution is used."
    ),
)
parser.add_argument(
    "--resolutionPrefitUncertaintyA",
    type=float,
    default=None,
    help=(
        "Override the relative prefit uncertainty for only the fitted "
        "resolution a parameter. If omitted, --resolutionPrefitUncertainty "
        "is used for a, b, c, and d."
    ),
)
parser = parsing.set_parser_default(parser, "theoryCorr", [])
parser = parsing.set_parser_default(parser, "scale_A", 5.0)
parser = parsing.set_parser_default(parser, "scale_e", 5.0)
parser = parsing.set_parser_default(parser, "fitMuonScaleAndResolution", True)

args = parser.parse_args()

if args.etaBins is not None and args.etaBins <= 0:
    raise ValueError("--etaBins must be a positive integer")
if args.resolutionPrefitUncertainty <= 0:
    raise ValueError("--resolutionPrefitUncertainty must be positive")
if (
    args.resolutionPrefitUncertaintyA is not None
    and args.resolutionPrefitUncertaintyA <= 0
):
    raise ValueError("--resolutionPrefitUncertaintyA must be positive")
if args.muonScaleVariation not in ["onnxReweight", "smearingWeightsSplines"]:
    raise ValueError(
        "d0_mass.py currently supports --muonScaleVariation onnxReweight "
        "or smearingWeightsSplines"
    )

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

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
    raise ValueError("Scale systematics require --muonCorrData massfit or lbl_massfit")

_smearing_helper, smearing_uncertainty_helper = (
    (None, None)
    if args.noSmearing
    else muon_calibration.make_muon_smearing_helpers(
        scale_var_method=args.muonScaleVariation,
        parameter_variations=True,
        fit_muon_resolution=args.fitMuonScaleAndResolution,
        variation_eta_bins=args.etaBins,
        resolution_prefit_uncertainties=[
            args.resolutionPrefitUncertaintyA or args.resolutionPrefitUncertainty,
            args.resolutionPrefitUncertainty,
            args.resolutionPrefitUncertainty,
            args.resolutionPrefitUncertainty,
        ],
        resolution_prefit_uncertainties_mode="relative",
    )
)

ROOT.gInterpreter.Declare("""
    #include <cmath>
    #include "ROOT/RVec.hxx"

    namespace wrem {
    namespace d0 {

    template <typename T>
    inline double first_or_value(const T &value) {
        return static_cast<double>(value);
    }

    template <typename T>
    inline double first_or_value(const ROOT::VecOps::RVec<T> &values) {
        return values.empty() ? -999.0 : static_cast<double>(values[0]);
    }

    }
    }
    """)


def limited_files(files):
    if args.maxFiles is not None and args.maxFiles > 0:
        return files[: args.maxFiles]
    return files


datasets = [
    narf.Dataset(
        name="D0_data",
        filepaths=limited_files([args.dataFile]),
        is_data=True,
        xsec=None,
        group="Data",
    ),
    narf.Dataset(
        name="D0_mc",
        filepaths=limited_files([args.mcFile]),
        is_data=False,
        xsec=1.0,
        group="D0",
    ),
]

axis_etaK = hist.axis.Regular(24, -2.4, 2.4, name="etaK")
axis_mRK = hist.axis.Variable(
    [
        0.00,
        0.06,
        0.10,
        0.14,
        0.18,
        0.22,
        0.26,
        0.30,
        0.34,
        0.38,
        0.42,
        0.46,
        0.50,
        0.60,
        0.70,
        0.85,
        1.05,
        1.30,
        1.55,
        1.80,
    ],
    name="mRK",
)
axis_D0mass = hist.axis.Regular(25, 1.8, 1.93, name="D0mass")
d0_axes = [axis_etaK, axis_mRK, axis_D0mass]

M_K = 0.493677
M_PI = 0.139570


def bool_filter(expression):
    return f"static_cast<bool>({expression})"


def scalarize(df, names):
    for name in names:
        df = df.Define(f"{name}0", f"wrem::d0::first_or_value({name})")
    return df


def define_reco(df):
    reco_cols = [
        "K_CVH_pt",
        "K_CVH_eta",
        "K_CVH_phi",
        "K_charge",
        "pi_CVH_pt",
        "pi_CVH_eta",
        "pi_CVH_phi",
        "pi_charge",
        "D0_CVH_mass",
    ]
    df = scalarize(df, reco_cols)
    return (
        df.Define(
            "K_E0",
            f"std::sqrt(std::pow(K_CVH_pt0*std::cosh(K_CVH_eta0), 2) + {M_K}*{M_K})",
        )
        .Define(
            "pi_E0",
            f"std::sqrt(std::pow(pi_CVH_pt0*std::cosh(pi_CVH_eta0), 2) + {M_PI}*{M_PI})",
        )
        .Define("D0_mass", "D0_CVH_mass0")
        .Define("mRK", f"{M_K}*{M_K}*(pi_E0/K_E0)/D0_mass")
    )


def reco_selection():
    return (
        "K_CVH_pt0 > 0.0 && pi_CVH_pt0 > 0.0 && K_E0 > 0.0 && pi_E0 > 0.0 && "
        "std::fabs(K_CVH_eta0) < 2.4 && "
        "D0_mass > 1.75 && D0_mass < 2.0"
    )


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")

    results = []
    df = df.DefinePerSample("weight", "1.0")
    weightsum = df.SumAndCount("weight")
    df = define_reco(df)
    df = df.Filter(bool_filter(reco_selection()))

    if dataset.is_data:
        results.append(
            df.HistoBoost(
                "hD0_data",
                d0_axes,
                ["K_CVH_eta0", "mRK", "D0_mass", "weight"],
            )
        )
    else:
        results.append(
            df.HistoBoost(
                "hD0_nom",
                d0_axes,
                ["K_CVH_eta0", "mRK", "D0_mass", "weight"],
            )
        )

        df = (
            df.Define("nominal_weight", "weight")
            .Define(
                "scale_recoPt",
                "ROOT::VecOps::RVec<float>{float(K_CVH_pt0), float(pi_CVH_pt0)}",
            )
            .Define(
                "scale_recoEta",
                "ROOT::VecOps::RVec<float>{float(K_CVH_eta0), float(pi_CVH_eta0)}",
            )
            .Define(
                "scale_recoPhi",
                "ROOT::VecOps::RVec<float>{float(K_CVH_phi0), float(pi_CVH_phi0)}",
            )
            .Define(
                "scale_recoCharge",
                "ROOT::VecOps::RVec<int>{int(K_charge0), int(pi_charge0)}",
            )
            # D0 ntuples do not carry gen-matched track kinematics; use the
            # reconstructed K/pi legs as the response-map coordinates.
            .Define(
                "scale_genPt",
                "ROOT::VecOps::RVec<float>{float(K_CVH_pt0), float(pi_CVH_pt0)}",
            )
            .Define(
                "scale_genEta",
                "ROOT::VecOps::RVec<float>{float(K_CVH_eta0), float(pi_CVH_eta0)}",
            )
            .Define(
                "scale_genPhi",
                "ROOT::VecOps::RVec<float>{float(K_CVH_phi0), float(pi_CVH_phi0)}",
            )
            .Define(
                "scale_genCharge",
                "ROOT::VecOps::RVec<int>{int(K_charge0), int(pi_charge0)}",
            )
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
                d0_axes,
                [
                    "K_CVH_eta0",
                    "mRK",
                    "D0_mass",
                    "nominal_muonScaleSyst_responseWeights_tensor",
                ],
                tensor_axes=data_jpsi_crctn_unc_helper.tensor_axes,
                storage=hist.storage.Double(),
            )
        )

        df = muon_calibration.add_resolution_uncertainty(
            df,
            d0_axes,
            results,
            ["K_CVH_eta0", "mRK", "D0_mass"],
            smearing_uncertainty_helper,
            "scale",
            storage_type=hist.storage.Double(),
            response_weight_col="scale_resolution_response_weight",
        )

    return results, weightsum


resultdict = narf.build_and_run(datasets, build_graph, event_tree=args.treeName)
fout = f"{os.path.basename(__file__).replace('py', 'hdf5')}"
write_analysis_output(resultdict, fout, args, name_append=[args.era])
