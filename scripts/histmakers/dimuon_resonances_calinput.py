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
    "--mcTriggerCut",
    default="1.",
    help="Additional cuts to apply to MC"
)
parser.add_argument(
    "--dataTriggerCut",
    default="1.",
    help="Additional cuts to apply to data"
)
parser = parsing.set_parser_default(parser, "theoryCorr", [])

args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

if args.muonScaleVariation not in ["onnxReweight", "smearingWeightsSplines"]:
    raise ValueError(
        "dimuon_resonances_calinput.py currently supports "
        "--muonScaleVariation onnxReweight or smearingWeightsSplines"
    )

calib_filepaths = common.calib_filepaths
diff_weights_helper = (
    ROOT.wrem.SplinesDifferentialWeightsHelper(calib_filepaths["tflite_file"])
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
    )
)


local_resonance_files = {
    "jpsi": {
        "data": ["/scratch/submit/cms/emanca/jpsicor_data.root"],
        "mc": [
            "/scratch/submit/cms/emanca/jpsicor_mc.root",
            "/scratch/submit/cms/emanca/jcor_mc_0to8.root",
        ],
    },
    # Temporary placeholder until the dedicated Upsilon ntuples are available.
    "upsilon": {
        "data": ["/scratch/submit/cms/emanca/jpsicor_data.root"],
        "mc": [
            "/scratch/submit/cms/emanca/jpsicor_mc.root",
            "/scratch/submit/cms/emanca/jcor_mc_0to8.root",
        ],
    },
}


def limited_files(files):
    if args.maxFiles is not None and args.maxFiles > 0:
        return files[: args.maxFiles]
    return files


datasets = [
    narf.Dataset(
        name=f"{args.resonance}_data",
        filepaths=limited_files(local_resonance_files[args.resonance]["data"]),
        is_data=True,
        group="Data",
    ),
    narf.Dataset(
        name=f"{args.resonance}_mc",
        filepaths=limited_files(local_resonance_files[args.resonance]["mc"]),
        is_data=False,
        xsec=1.0,
        group=args.resonance.upper(),
    ),
]


resonance_options = {
    "jpsi": {
        "name": "JPsi",
        "eta_axis": hist.axis.Regular(24, -2.4, 2.4, name="eta1"),
        "eta2_axis": hist.axis.Regular(24, -2.4, 2.4, name="eta2"),
        "pt_axis": hist.axis.Variable([4.2, 7.0, 10.5, 15.0, 25.0], name="pt1"),
        "pt2_axis": hist.axis.Variable([4.2, 7.0, 10.5, 15.0, 25.0], name="pt2"),
        "mass_axis": hist.axis.Regular(25, 2.92, 3.28, name="mass"),
        "selection": (
            "Mupluscor_pt > 1.0 && Muminuscor_pt > 1.0 && "
            "Mupluscor_pt < 100.0 && Muminuscor_pt < 100.0 && "
            "std::fabs(Mupluscor_eta) < 2.4 && std::fabs(Muminuscor_eta) < 2.4 && "
            "Jpsicor_mass > 2.8 && Jpsicor_mass < 3.35"
        ),
    },
    "upsilon": {
        "name": "Y",
        "eta_axis": hist.axis.Regular(8, -0.8, 0.8, name="eta1"),
        "eta2_axis": hist.axis.Regular(8, -0.8, 0.8, name="eta2"),
        "pt_axis": hist.axis.Variable([4.2, 6.0, 7.9, 10.3, 25.0], name="pt1"),
        "pt2_axis": hist.axis.Variable([4.2, 6.0, 7.9, 10.3, 25.0], name="pt2"),
        "mass_axis": hist.axis.Regular(25, 9.0, 9.7, name="mass"),
        "selection": (
            "Mupluscor_pt > 1.0 && Muminuscor_pt > 1.0 && "
            "Mupluscor_pt < 100.0 && Muminuscor_pt < 100.0 && "
            "std::fabs(Mupluscor_eta) < 0.8 && std::fabs(Muminuscor_eta) < 0.8 && "
            "Jpsicor_mass > 8.8 && Jpsicor_mass < 9.6"
        ),
    },
}

cfg = resonance_options[args.resonance]
calibration_axes = [
    cfg["eta_axis"],
    cfg["eta2_axis"],
    cfg["pt_axis"],
    cfg["pt2_axis"],
    cfg["mass_axis"],
]
calibration_cols = [
    "Mupluscor_eta",
    "Muminuscor_eta",
    "Mupluscor_pt",
    "Muminuscor_pt",
    "Jpsicor_mass",
]


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")

    results = []

    df = df.DefinePerSample("weight", "1.0")
    weightsum = df.SumAndCount("weight")

    df = df.Filter(cfg["selection"])
    if not dataset.is_data:
        df = df.Filter(
            "Jpsigen_mass > 0.0 && Muplusgen_pt > 0.0 && Muminusgen_pt > 0.0"
        )
        df = df.Filter(
            args.mcTriggerCut
        )
    else:
        df = df.Filter(
            args.dataTriggerCut
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
                "ROOT::VecOps::RVec<float>{float(Mupluscor_pt), float(Muminuscor_pt)}",
            )
            .Define(
                "scale_recoEta",
                "ROOT::VecOps::RVec<float>{float(Mupluscor_eta), float(Muminuscor_eta)}",
            )
            .Define(
                "scale_recoPhi",
                "ROOT::VecOps::RVec<float>{float(Mupluscor_phi), float(Muminuscor_phi)}",
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
        if diff_weights_helper is not None:
            df = df.Define(
                "scale_response_weight",
                diff_weights_helper,
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
        )

    return results, weightsum


resultdict = narf.build_and_run(datasets, build_graph, event_tree="tree")

fout = f"{os.path.basename(__file__).replace('py', 'hdf5')}"
write_analysis_output(resultdict, fout, args, name_append=[args.resonance, args.era])
