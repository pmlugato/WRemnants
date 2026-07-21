import glob
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
    default="/scratch/submit/cms/emanca/D0Data_v4.root",
    help="Input ROOT file for data",
)
parser.add_argument(
    "--mcFile",
    default="/scratch/submit/cms/emanca/D0MC",
    help="Input ROOT file, glob, or directory for MC",
)
parser.add_argument(
    "--treeName",
    default="tree",
    help="Input tree name",
)
parser.add_argument(
    "--applyFitDeltaM",
    action="store_true",
    help=(
        "Also apply the later fitted DeltaM cut "
        "abs(deltaM_D0pis_piK - 0.14543) < 0.003."
    ),
)
parser.add_argument(
    "--pvalCut",
    type=float,
    default=0.005,
    help=(
        "Minimum vertex-fit p-value required for both the piK (D0) and D0pis (D*) "
        "vertices. Disabled by default: the current ntuples store an unfilled -99 "
        "sentinel in pval_piK/pval_D0pis, so any positive cut would reject every "
        "event. Set e.g. 0.005 once valid p-values are produced. When enabled the "
        "cut is guarded so sentinel (-99) values are never rejected, and the "
        "pval_piK/pval_D0pis branches are only required in that case."
    ),
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
    default=0.3,
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
parser.add_argument(
    "--buildManualScaleVariations",
    action="store_true",
    help="Include manual scale variations in addition to reweights"
)
parser = parsing.set_parser_default(parser, "theoryCorr", [])
parser = parsing.set_parser_default(parser, "scale_A", 1.0)
parser = parsing.set_parser_default(parser, "scale_e", 1.0)
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
    #include <algorithm>
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

    inline double delta_phi(double phi1, double phi2) {
        double dphi = phi1 - phi2;
        while (dphi > M_PI) dphi -= 2.0*M_PI;
        while (dphi <= -M_PI) dphi += 2.0*M_PI;
        return dphi;
    }

    template <typename PtVec, typename EtaVec, typename PhiVec, typename ChargeVec, typename PdgIdVec>
    ROOT::VecOps::RVec<float> matched_gen_kinematics(
        double reco_pt,
        double reco_eta,
        double reco_phi,
        int reco_charge,
        const PtVec &gen_pt,
        const EtaVec &gen_eta,
        const PhiVec &gen_phi,
        const ChargeVec &gen_charge,
        const PdgIdVec &gen_pdgId,
        int abs_pdgid,
        double max_dr = 0.03
    ) {
        const auto size = std::min({gen_pt.size(), gen_eta.size(), gen_phi.size(),
                                    gen_charge.size(), gen_pdgId.size()});
        const double max_dr2 = max_dr*max_dr;
        double best_dr2 = max_dr2;
        int best = -1;
        for (std::size_t i = 0; i < size; ++i) {
            if (std::abs(static_cast<int>(gen_pdgId[i])) != abs_pdgid) continue;
            if (static_cast<int>(gen_charge[i]) != reco_charge) continue;

            const double deta = reco_eta - static_cast<double>(gen_eta[i]);
            const double dphi = delta_phi(reco_phi, static_cast<double>(gen_phi[i]));
            const double dr2 = deta*deta + dphi*dphi;
            if (dr2 < best_dr2) {
                best_dr2 = dr2;
                best = static_cast<int>(i);
            }
        }

        if (best < 0) return ROOT::VecOps::RVec<float>{0.f, 0.f, 0.f, 0.f, 0.f};

        return ROOT::VecOps::RVec<float>{
            static_cast<float>(gen_pt[best]),
            static_cast<float>(gen_eta[best]),
            static_cast<float>(gen_phi[best]),
            static_cast<float>(gen_charge[best]),
            1.f
        };
    }

    }
    }
    """)

# Manual A/e/M scale variations: recompute the D0 mass (and mRK) from the K/pi
# four-vectors under an independent +-1 prefit-width shift of each (eta bin, scale
# parameter), instead of reweighting. The correction model matches
# muon_calibration.define_AeM_data_corrections:
#     pt_cor = (1 + A' - e'/pt + q*M'*pt) * pt   (primed = physical units)
# The D0 mass is linear in each of A', e', M', so the per-variation pt shift is
# independent of the (unknown here) central A/e/M values and we propagate the shift
# directly:  dpt_A = dA'*pt,  dpt_e = -de',  dpt_M = q*dM'*pt^2.
# Both the varied mass and mRK are expressed as a shift on top of their authoritative
# nominal values (D0_CVH_mass0 and the nominal mRK column), so a variation that touches
# neither track's eta bin lands exactly on the nominal histogram. Mass and mRK share the
# same shifted momenta (computed once) so they never drift apart.
ROOT.gInterpreter.Declare("""
    #include <cmath>
    #include "ROOT/RVec.hxx"
    #include "Math/Vector4D.h"

    namespace wrem {
    namespace d0 {

    // Per-event varied observables: mass[k] and mRK[k] share the flattened index k.
    struct ScaleVarResult {
        ROOT::VecOps::RVec<double> mass;
        ROOT::VecOps::RVec<double> mRK;
    };

    class ScaleVariationsHelper {
    public:
        ScaleVariationsHelper(unsigned int n_eta_bins, double eta_min, double eta_max,
                              unsigned int n_scale_params,
                              double width_A, double width_e, double width_M,
                              double mK, double mPi)
          : n_eta_bins_(n_eta_bins), n_scale_params_(n_scale_params),
            nvars_(n_eta_bins * n_scale_params),
            eta_min_(eta_min), eta_max_(eta_max),
            width_A_(width_A), width_e_(width_e), width_M_(width_M),
            mK_(mK), mPi_(mPi) {}

        int eta_bin(double eta) const {
            if (eta < eta_min_ || eta >= eta_max_) return -1;
            int ib = static_cast<int>((eta - eta_min_) / (eta_max_ - eta_min_) * n_eta_bins_);
            if (ib < 0 || ib >= static_cast<int>(n_eta_bins_)) return -1;
            return ib;
        }

        // pt of a track after shifting scale parameter iparm in eta bin ieta by sgn*width
        double shifted_pt(double pt, double q, int track_ieta,
                          unsigned int ieta, unsigned int iparm, double sgn) const {
            if (track_ieta != static_cast<int>(ieta)) return pt;
            double dpt = 0.0;
            if (iparm == 0)      dpt = sgn * width_A_ * pt;            // A': dpt = dA'*pt
            else if (iparm == 1) dpt = -sgn * width_e_;               // e': dpt = -de'
            else if (iparm == 2) dpt = sgn * width_M_ * q * pt * pt;  // M': dpt = q*dM'*pt^2
            return pt + dpt;
        }

        // mRK observable, matching the nominal column M_K^2 * (pi_E/K_E) / mass.
        double mrk(double eK, double ePi, double mass) const {
            return mK_ * mK_ * (ePi / eK) / mass;
        }

        // Returns nvars*2 varied mass and mRK values, ordered (ieta, iparm, {down, up})
        // so the flattened index is ivar*2 + isign with ivar = ieta*n_scale_params + iparm.
        ScaleVarResult operator()(
            double ptK, double etaK, double phiK, double qK,
            double ptPi, double etaPi, double phiPi, double qPi,
            double d0_mass_nom, double mRK_nom) const {

            const int biK = eta_bin(etaK);
            const int biPi = eta_bin(etaPi);

            const double coshK = std::cosh(etaK);
            const double coshPi = std::cosh(etaPi);
            const auto energy = [](double pt, double cosh_eta, double m) {
                return std::sqrt(pt * cosh_eta * pt * cosh_eta + m * m);
            };

            const ROOT::Math::PtEtaPhiMVector v0K(ptK, etaK, phiK, mK_);
            const ROOT::Math::PtEtaPhiMVector v0Pi(ptPi, etaPi, phiPi, mPi_);
            const double m_nom = (v0K + v0Pi).M();
            const double mrk_nom_recomp =
                mrk(energy(ptK, coshK, mK_), energy(ptPi, coshPi, mPi_), m_nom);

            ScaleVarResult res;
            res.mass.reserve(nvars_ * 2);
            res.mRK.reserve(nvars_ * 2);
            const double signs[2] = {-1.0, 1.0};  // down, up
            for (unsigned int ieta = 0; ieta < n_eta_bins_; ++ieta) {
                for (unsigned int iparm = 0; iparm < n_scale_params_; ++iparm) {
                    for (int isign = 0; isign < 2; ++isign) {
                        const double sgn = signs[isign];
                        const double ptKs = shifted_pt(ptK, qK, biK, ieta, iparm, sgn);
                        const double ptPis = shifted_pt(ptPi, qPi, biPi, ieta, iparm, sgn);
                        const ROOT::Math::PtEtaPhiMVector vK(ptKs, etaK, phiK, mK_);
                        const ROOT::Math::PtEtaPhiMVector vPi(ptPis, etaPi, phiPi, mPi_);
                        const double m_var = (vK + vPi).M();
                        const double mrk_var =
                            mrk(energy(ptKs, coshK, mK_), energy(ptPis, coshPi, mPi_), m_var);
                        res.mass.push_back(d0_mass_nom + (m_var - m_nom));
                        res.mRK.push_back(mRK_nom + (mrk_var - mrk_nom_recomp));
                    }
                }
            }
            return res;
        }

    private:
        unsigned int n_eta_bins_, n_scale_params_, nvars_;
        double eta_min_, eta_max_, width_A_, width_e_, width_M_, mK_, mPi_;
    };

    // Constant per-fill index vectors matching the (ivar, {down, up}) ordering above.
    inline ROOT::VecOps::RVec<int> scale_var_indices(unsigned int nvars) {
        ROOT::VecOps::RVec<int> out;
        out.reserve(nvars * 2);
        for (unsigned int ivar = 0; ivar < nvars; ++ivar) {
            out.push_back(static_cast<int>(ivar));
            out.push_back(static_cast<int>(ivar));
        }
        return out;
    }
    // down -> -1.0 (downUpVar bin 0), up -> +1.0 (downUpVar bin 1)
    inline ROOT::VecOps::RVec<double> updown_coords(unsigned int nvars) {
        ROOT::VecOps::RVec<double> out;
        out.reserve(nvars * 2);
        for (unsigned int ivar = 0; ivar < nvars; ++ivar) {
            out.push_back(-1.0);
            out.push_back(1.0);
        }
        return out;
    }

    }
    }
    """)


def limited_files(files):
    if args.maxFiles is not None and args.maxFiles > 0:
        return files[: args.maxFiles]
    return files


def input_files(path):
    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, "**", "*.root"), recursive=True))
    else:
        files = sorted(glob.glob(path))
    if not files:
        raise ValueError(f"No input ROOT files found for {path}")
    return limited_files(files)


datasets = [
    narf.Dataset(
        name="D0_data",
        filepaths=input_files(args.dataFile),
        is_data=True,
        xsec=None,
        group="Data",
    ),
    narf.Dataset(
        name="D0_mc",
        filepaths=input_files(args.mcFile),
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

# Manual A/e/M scale-variation histogram. The per-(eta bin, parameter) construction
# below is only valid for the diagonal prefit-width configuration (no covariance /
# eigen-decomposition), i.e. --fitMuonScaleAndResolution, which is the default here.
# The scale-variation and down/up axes are reused verbatim from the response-weight
# uncertainty helper so the manual histogram is structurally identical to
# nominal_muonScaleSyst_responseWeights and the "unc" index maps to the same nuisance.
build_manual_scale_variations = args.buildManualScaleVariations
d0_scale_var_helper = None
manual_scale_axes = None
if build_manual_scale_variations:
    _scale_var_axis, _updown_axis = data_jpsi_crctn_unc_helper.tensor_axes
    _nvars = _scale_var_axis.size
    if _nvars % 3 != 0:
        raise ValueError(
            f"Expected the scale unc axis size to be a multiple of 3 (A, e, M per "
            f"eta bin), got {_nvars}"
        )
    _n_eta_bins_scale = _nvars // 3
    d0_scale_var_helper = ROOT.wrem.d0.ScaleVariationsHelper(
        _n_eta_bins_scale,
        -2.4,
        2.4,
        3,
        1e-3 * args.scale_A,
        1e-2 * args.scale_e,
        1e-4 * args.scale_M,
        M_K,
        M_PI,
    )
    manual_scale_axes = [_scale_var_axis, _updown_axis]
    logger.info(
        f"Manual A/e/M scale variations enabled: {_nvars} nuisances "
        f"({_n_eta_bins_scale} eta bins x 3 parameters), down/up each."
    )
else:
    logger.warning(
        "Manual A/e/M scale variations require --fitMuonScaleAndResolution "
        "(diagonal prefit widths); skipping nominal_muonScaleSyst_manual."
    )


def bool_filter(expression):
    return f"static_cast<bool>({expression})"


def scalarize(df, names):
    for name in names:
        df = df.Define(f"{name}0", f"wrem::d0::first_or_value({name})")
    return df


def available_columns(df):
    return {str(col) for col in df.GetColumnNames()}


def require_columns(columns, names):
    missing = [name for name in names if name not in columns]
    if missing:
        raise ValueError(f"Missing required D0 ntuple branches: {', '.join(missing)}")


def available_first(columns, names, label):
    available = [name for name in names if name in columns]
    if not available:
        raise ValueError(
            f"Missing D0 ntuple branch for {label}; tried {', '.join(names)}"
        )
    return available


def define_reco(df):
    columns = available_columns(df)
    require_columns(
        columns,
        [
            "K_CVH_pt",
            "K_CVH_eta",
            "K_CVH_phi",
            "K_charge",
            "pi_CVH_pt",
            "pi_CVH_eta",
            "pi_CVH_phi",
            "pi_charge",
        ],
    )
    # The vertex-fit p-value branches are only needed when the cut is enabled;
    # they are unfilled (-99) in current productions, so do not require them by default.
    pval_cols = ["pval_piK", "pval_D0pis"] if args.pvalCut is not None else []
    if pval_cols:
        require_columns(columns, pval_cols)
    d0_mass_cols = available_first(
        columns, ["mass_piK", "D0_fit_mass", "D0_CVH_mass"], "D0 fitted mass"
    )
    dst_mass_cols = available_first(
        columns,
        ["mass_D0pis", "Dst_fit_mass", "Dst_CVH3_mass", "Dst_CVH_mass"],
        "D* fitted mass",
    )
    delta_m_cols = [
        name for name in ["deltaM_D0pis_piK", "Dst_fit_deltaM"] if name in columns
    ]

    reco_cols = [
        "K_CVH_pt",
        "K_CVH_eta",
        "K_CVH_phi",
        "K_charge",
        "pi_CVH_pt",
        "pi_CVH_eta",
        "pi_CVH_phi",
        "pi_charge",
        *pval_cols,
        *d0_mass_cols,
        *dst_mass_cols,
        *delta_m_cols,
    ]
    df = scalarize(df, reco_cols)
    d0_expr = f"{d0_mass_cols[-1]}0"
    for name in reversed(d0_mass_cols[:-1]):
        d0_expr = f"({name}0 > 0.0 ? {name}0 : {d0_expr})"
    dst_expr = f"{dst_mass_cols[-1]}0"
    for name in reversed(dst_mass_cols[:-1]):
        dst_expr = f"({name}0 > 0.0 ? {name}0 : {dst_expr})"
    if delta_m_cols:
        delta_m_expr = f"{delta_m_cols[-1]}0"
        for name in reversed(delta_m_cols[:-1]):
            delta_m_expr = f"({name}0 > 0.0 ? {name}0 : {delta_m_expr})"
    else:
        delta_m_expr = "Dst_fit_mass_for_selection - D0_fit_mass_for_selection"

    df = (
        df.Define(
            "K_E0",
            f"std::sqrt(std::pow(K_CVH_pt0*std::cosh(K_CVH_eta0), 2) + {M_K}*{M_K})",
        )
        .Define(
            "pi_E0",
            f"std::sqrt(std::pow(pi_CVH_pt0*std::cosh(pi_CVH_eta0), 2) + {M_PI}*{M_PI})",
        )
        .Define("D0_fit_mass_for_selection", d0_expr)
        .Define("Dst_fit_mass_for_selection", dst_expr)
        .Define("Dst_fit_deltaM_for_selection", delta_m_expr)
        .Define("D0_mass", "D0_fit_mass_for_selection")
        .Define("mRK", f"{M_K}*{M_K}*(pi_E0/K_E0)/D0_mass")
    )
    if args.applyFitDeltaM and not delta_m_cols:
        logger.warning(
            "No fitted DeltaM branch found; using D* fitted mass minus D0 fitted mass"
        )
    return df


def reco_selection():
    selection = [
        "K_CVH_pt0 > 0.0",
        "pi_CVH_pt0 > 0.0",
        "K_E0 > 0.0",
        "pi_E0 > 0.0",
        "std::fabs(K_CVH_eta0) < 2.4",
        "D0_fit_mass_for_selection > 0.0",
        "Dst_fit_mass_for_selection > 0.0",
        "std::fabs(D0_fit_mass_for_selection - 1.86483) < 0.05",
        "std::fabs(Dst_fit_mass_for_selection - 2.01026) < 0.15",
    ]
    if args.pvalCut is not None:
        # Guard against the unfilled -99 sentinel: reject a candidate only when its
        # p-value is valid (> -1) and below threshold. Sentinel values are kept.
        selection.append(f"(pval_piK0 < -1.0 || pval_piK0 > {args.pvalCut})")
        selection.append(f"(pval_D0pis0 < -1.0 || pval_D0pis0 > {args.pvalCut})")
    if args.applyFitDeltaM:
        selection.append("std::fabs(Dst_fit_deltaM_for_selection - 0.14543) < 0.003")
    return " && ".join(selection)


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
            .Define(
                "K_gen_match",
                "wrem::d0::matched_gen_kinematics("
                "K_CVH_pt0, K_CVH_eta0, K_CVH_phi0, int(K_charge0), "
                "gen_pt, gen_eta, gen_phi, gen_charge, gen_pdgId, 321)",
            )
            .Define(
                "pi_gen_match",
                "wrem::d0::matched_gen_kinematics("
                "pi_CVH_pt0, pi_CVH_eta0, pi_CVH_phi0, int(pi_charge0), "
                "gen_pt, gen_eta, gen_phi, gen_charge, gen_pdgId, 211)",
            )
            .Define(
                "scale_genPt",
                "ROOT::VecOps::RVec<float>{K_gen_match[0], pi_gen_match[0]}",
            )
            .Define(
                "scale_genEta",
                "ROOT::VecOps::RVec<float>{K_gen_match[1], pi_gen_match[1]}",
            )
            .Define(
                "scale_genPhi",
                "ROOT::VecOps::RVec<float>{K_gen_match[2], pi_gen_match[2]}",
            )
            .Define(
                "scale_genCharge",
                "ROOT::VecOps::RVec<int>{int(K_gen_match[3]), int(pi_gen_match[3])}",
            )
            .Define(
                "scale_hasGenMatch",
                "ROOT::VecOps::RVec<int>{int(K_gen_match[4]), int(pi_gen_match[4])}",
            )
            .Define("scale_muon_source", "ROOT::VecOps::RVec<int>{443, 443}")
            .Filter("K_gen_match[4] > 0.5 && pi_gen_match[4] > 0.5")
        )

        results.append(
            df.HistoBoost(
                "hD0_nom",
                d0_axes,
                ["K_CVH_eta0", "mRK", "D0_mass", "weight"],
            )
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

        if build_manual_scale_variations:
            # Manual scale variations: recompute the D0 mass and mRK under each +-1
            # prefit-width A/e/M shift and fill a single histogram in one pass. The mass
            # and mRK coordinates are RVecs of length nvars*2 (both varied, sharing the
            # same shifted momenta); etaK is held at its nominal reco value since a pt
            # scale does not change eta.
            df = df.Define(
                "d0_scale_var",
                d0_scale_var_helper,
                [
                    "K_CVH_pt0",
                    "K_CVH_eta0",
                    "K_CVH_phi0",
                    "K_charge0",
                    "pi_CVH_pt0",
                    "pi_CVH_eta0",
                    "pi_CVH_phi0",
                    "pi_charge0",
                    "D0_mass",
                    "mRK",
                ],
            )
            df = df.Define("d0_scale_var_mass", "d0_scale_var.mass")
            df = df.Define("d0_scale_var_mRK", "d0_scale_var.mRK")
            df = df.Define(
                "d0_scale_var_idx",
                f"wrem::d0::scale_var_indices({manual_scale_axes[0].size})",
            )
            df = df.Define(
                "d0_scale_updown",
                f"wrem::d0::updown_coords({manual_scale_axes[0].size})",
            )
            results.append(
                df.HistoBoost(
                    "nominal_muonScaleSyst_manual",
                    d0_axes + manual_scale_axes,
                    [
                        "K_CVH_eta0",
                        "d0_scale_var_mRK",
                        "d0_scale_var_mass",
                        "d0_scale_var_idx",
                        "d0_scale_updown",
                        "nominal_weight",
                    ],
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
