import os

import hist
import ROOT

import narf
from wremnants.production import muon_calibration
from wremnants.production.histmaker_tools import write_analysis_output
from wremnants.utilities import common, parsing
from wums import logging

analysis_label = common.analysis_label(os.path.basename(__file__))


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


parser, initargs = parsing.common_parser(analysis_label)
parser.add_argument(
    "--resonance",
    type=str,
    choices=["jpsi", "upsilon"],
    help="Resonance for selection",
    default="upsilon",
)
parser.add_argument(
    "--triggers",
    nargs="+",
    default=None,
    help="Triggers to process for the resonance (each as a separate channel)",
)
parser.add_argument(
    "--etaBins",
    type=int,
    default=None,
    help="Override the number of bins for each reconstructed muon eta axis",
)
parser.add_argument(
    "--ptBinEdges",
    type=float,
    nargs="+",
    default=None,
    help="Override the default pt bining for the selected resonance"
)
parser.add_argument(
    "--applyAeMtoData",
    type=str,
    default=None,
    help="Path to a rabbit fitresult HDF5 file. The fitted AeM parameters are "
    "extracted and used to correct the data muon pts "
    "(via muon_calibration.define_AeM_data_corrections) before filling the "
    "data histograms. MC histograms are left unchanged.",
)
parser = parsing.set_parser_default(parser, "theoryCorr", [])

args = parser.parse_args()

if args.etaBins is not None and args.etaBins <= 0:
    raise ValueError("--etaBins must be a positive integer")
if args.etaBins is not None and not args.fitMuonScaleAndResolution:
    raise ValueError(
        "--etaBins currently requires --fitMuonScaleAndResolution so scale "
        "uncertainties can be correlated across the requested eta groups"
    )
possible_triggers = [trig["label"] for trig in trigger_channels[args.resonance]]
selected_triggers = possible_triggers if args.triggers is None else args.triggers
if not all(trig in possible_triggers for trig in selected_triggers):
    raise ValueError(
        f"--triggers must contain only valid triggers for the "
        f"selected resonance ({possible_triggers} for resonance "
        f"{args.resonance}), or be None to use all valid triggers"
    )

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)


def read_AeM_corrections(fitresult_path, n_eta_bins=24):
    import rabbit.io_tools

    fitresult = rabbit.io_tools.get_fitresult(fitresult_path)
    h_parms = fitresult["parms"].get()
    labels = [str(label) for label in h_parms.axes["parms"]]
    values = h_parms.values()
    parm_map = dict(zip(labels, values))

    def collect(prefix):
        out = []
        for i in range(n_eta_bins):
            key = f"{prefix}{i}"
            if key not in parm_map:
                raise ValueError(
                    f"Parameter '{key}' not found in fitresult '{fitresult_path}'. "
                    f"Available parameters: {labels}"
                )
            out.append(float(parm_map[key]))
        return out

    return collect("A"), collect("e"), collect("M")


AeM_corrections = None
if args.applyAeMtoData is not None:
    AeM_corrections = read_AeM_corrections(args.applyAeMtoData)
    logger.info(
        f"Loaded A/e/M data pt corrections from fitresult {args.applyAeMtoData}"
    )

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
        resolution_prefit_uncertainties=[0.3, 0.3, 0.3, 0.3],
        resolution_prefit_uncertainties_mode="relative",
    )
)

(
    _pixel_multiplicity_helper,
    pixel_multiplicity_uncertainty_helper,
    pixel_multiplicity_uncertainty_helper_stat,
) = muon_calibration.make_pixel_multiplicity_helpers()


def limited_files(files):
    if args.maxFiles is not None and args.maxFiles > 0:
        return files[: args.maxFiles]
    return files


def bool_filter(expression):
    return f"static_cast<bool>({expression})"


datasets = []
dataset_channels = {}
for channel in trigger_channels[args.resonance]:
    if channel["label"] not in selected_triggers:
        continue
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
        "pt1_axis": hist.axis.Variable([4.2, 7.0, 10.5, 15.0, 25.0], name="pt1"),
        "pt2_axis": hist.axis.Variable([4.2, 7.0, 10.5, 15.0, 25.0], name="pt2"),
        "mass_axis": hist.axis.Regular(25, 2.92, 3.28, name="mass"),
        "mass_range": (2.8, 3.35),
    },
    "upsilon": {
        "name": "Y",
        "default_eta_bins": 8,
        "eta_range": (-0.8, 0.8),
        "pt1_axis": hist.axis.Variable([4.2, 6.0, 7.9, 10.3, 25.0], name="pt1"),
        "pt2_axis": hist.axis.Variable([4.2, 6.0, 7.9, 10.3, 25.0], name="pt2"),
        "mass_axis": hist.axis.Regular(25, 9.0, 9.7, name="mass"),
        "mass_range": (8.8, 9.6),
    },
}

cfg = resonance_options[args.resonance]
eta_bins = args.etaBins or cfg["default_eta_bins"]
eta_min, eta_max = cfg["eta_range"]
calibration_axes = [
    hist.axis.Regular(eta_bins, eta_min, eta_max, name="eta1"),
    hist.axis.Regular(eta_bins, eta_min, eta_max, name="eta2"),
    cfg["pt1_axis"] if args.ptBinEdges is None else hist.axis.Variable(args.ptBinEdges, name="pt1"),
    cfg["pt2_axis"] if args.ptBinEdges is None else hist.axis.Variable(args.ptBinEdges, name="pt2"),
    cfg["mass_axis"],
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


def build_graph(df, dataset):
    logger.info(f"build graph for dataset: {dataset.name}")

    results = []
    channel = dataset_channels[dataset.name]
    reco_cols = reco_columns(channel)

    if dataset.is_data and AeM_corrections is not None:
        if not channel["layer_corrected"]:
            raise ValueError(
                "--applyAeMtoData operates on the layer-corrected data columns "
                f"(Mupluscor_*), but data channel '{channel['label']}' is not "
                "layer corrected"
            )
        A_vals, e_vals, M_vals = AeM_corrections
        df = muon_calibration.define_AeM_data_corrections(df, A_vals, e_vals, M_vals)
        reco_cols["plus_pt"] = "Mupluscor_AeM_pt"
        reco_cols["minus_pt"] = "Muminuscor_AeM_pt"
        reco_cols["mass"] = "Jpsicor_AeM_mass"

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
    if not dataset.is_data:
        df = df.Filter(
            "Jpsigen_mass > 0.0 && Muplusgen_pt > 0.0 && Muminusgen_pt > 0.0"
        )
    df = df.Filter(bool_filter(channel["cut"]))

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

        pixel_multiplicity_cols = [
            "pixel_triggerCat",
            "pixel_eta",
            "pixel_pt",
            "pixel_charge",
            "pixel_nvalidpixel",
        ]

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
write_analysis_output(resultdict, fout, args, name_append=name_append)
