import argparse
import sys
from pathlib import Path

import h5py
import hist
import numpy as np

from rabbit import tensorwriter

WREMNANTS_DIR = Path(__file__).resolve().parent / "WRemnants"
if str(WREMNANTS_DIR) not in sys.path:
    sys.path.insert(0, str(WREMNANTS_DIR))

from wremnants.utilities.io_tools import base_io
from wums import ioutils

jpsi_channel_names = {
    "dimuon20_jpsi": "JPsi_prompt",
    "doublemu4_jpsitrk_displaced": "JPsi_displaced",
}

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", default="./", help="output directory")
parser.add_argument("--outname", default="test_tensor", help="output file name")
parser.add_argument(
    "--postfix",
    default=None,
    type=str,
    help="Postfix to append on output file name",
)
parser.add_argument(
    "--sparse",
    default=False,
    action="store_true",
    help="Make sparse tensor",
)
parser.add_argument(
    "--symmetrizeAll",
    default=False,
    action="store_true",
    help="Make fully symmetric tensor",
)
parser.add_argument(
    "--skipMaskedChannels",
    default=False,
    action="store_true",
    help="Skip adding masked channels",
)
parser.add_argument(
    "--systematicType",
    choices=["log_normal", "normal"],
    default="log_normal",
    help="probability density for systematic variations",
)
parser.add_argument(
    "--calinputHdf5",
    default=str(WREMNANTS_DIR / "dimuon_resonances_calinput_jpsi_2016PostVFP.hdf5"),
    help="WRemnants calInput-style HDF5 file to use for the JPsi channel",
)
parser.add_argument(
    "--jpsiChannels",
    nargs="+",
    choices=list(jpsi_channel_names.values()),
    default=list(jpsi_channel_names.values()),
    help="JPsi channels to add to the tensor",
)
parser.add_argument(
    "--mcEventThresh",
    default=-1,
    type=float,
    help="event threshold in raw MC for a 4d bin being included in the fit. Disabled by default. Zeros both MC and data in failing bins if specified",
)
parser.add_argument(
    "--clip",
    default=0,
    type=float,
    help="clip values in MC and data at this threshold",
)
args = parser.parse_args()


# rebin differently depending on which sample
# temperorarily lower mass bins bc numerical bs
def rebin(bh, sample):
    if sample == "JPsi":
        return bh[
            {
                "pt1": np.s_[::1j],
                "pt2": np.s_[::1j],
                "eta1": np.s_[::3j],
                "eta2": np.s_[::3j],
                "mass": np.s_[::1j],
            }
        ]
    raise Exception("invalid sample")


# modifies the input histogram variable and then returns it for convenience
def clip_values(dummy_hist, lim=1e-6, variances=True):
    # return dummy_hist #try getting rid of this
    values = np.copy(dummy_hist.values())
    values_copy = np.copy(np.clip(values, lim, np.inf))
    dummy_hist.values()[...] = values_copy
    if variances:
        variances = dummy_hist.variances()
        if variances is not None:
            variances_copy = np.copy(np.clip(variances, lim, np.inf))
            dummy_hist.variances()[...] = variances_copy
    return dummy_hist


def populated_cells(dummy_hist, shape_axis="mass", event_thresh=0):
    """Return the nonempty cells after integrating over the lineshape axis."""
    shape_axis_index = dummy_hist.axes.name.index(shape_axis)
    return (np.sum(dummy_hist.values(), axis=shape_axis_index) > max(event_thresh, 0))


def zero_inactive_cells(dummy_hist, active_cells, shape_axis="mass"):
    """Restore exact zeros outside active kinematic cells after clipping."""
    shape_axis_index = dummy_hist.axes.name.index(shape_axis)
    active_bins = np.expand_dims(active_cells, axis=shape_axis_index)
    dummy_hist.values()[...] = np.where(active_bins, dummy_hist.values(), 0.0)
    variances = dummy_hist.variances()
    if variances is not None:
        dummy_hist.variances()[...] = np.where(active_bins, variances, 0.0)
    return dummy_hist


def make_active_background(template_hist, active_cells, shape_axis="mass"):
    """Build a flat background only in populated kinematic cells."""
    shape_axis_index = template_hist.axes.name.index(shape_axis)
    active_bins = np.expand_dims(active_cells, axis=shape_axis_index)
    background_hist = hist.Hist(*template_hist.axes, storage=hist.storage.Weight())
    background_hist.values()[...] = np.broadcast_to(
        active_bins, template_hist.values().shape
    )
    background_hist.variances()[...] = 0
    return background_hist


# def hist_rescale(hist, lumi, xsec, mc_events):
#     return hist * lumi * xsec * 1000 / mc_events


# numerical stupidity can cause problems - fix by setting vars to nom if they differ by a rounding error
# modifies the input histogram and then returns it
def trim_variations(var, nom, thresh=1e-5):
    var_values = var.values()
    nom_values = nom.values()
    indices = np.where(np.abs(var_values - nom_values) < thresh)
    var.values()[indices] = nom_values[indices]
    return var


def materialize(obj):
    return obj.get() if isinstance(obj, ioutils.H5PickleProxy) else obj


def load_calinput_datasets(path):
    with h5py.File(path, "r") as h5file:
        results = base_io.load_results_h5py(h5file)

        datasets = {}
        for dataset_name, result in results.items():
            if dataset_name == "meta_info" or not isinstance(result, dict):
                continue
            histograms = {}
            for hist_name, hist_obj in result.get("output", {}).items():
                histograms[hist_name] = materialize(hist_obj)
            datasets[dataset_name] = histograms

    return datasets





def make_jpsi_channel(dataset_label, histograms):
    required = [
        "JPsi_data",
        "JPsi_mc",
        "nominal_muonScaleSyst_responseWeights",
        "nominal_muonResolutionSyst_responseWeights",
        "nominal_pixelMultiplicitySyst",
    ]
    missing = [name for name in required if name not in histograms]
    if missing:
        raise RuntimeError(
            f"Missing histograms for JPsi channel '{dataset_label}': {missing}. "
            "Run dimuon_resonances_calinput.py with MC scale/pixel systematics enabled."
        )

    channel = {
        "name": jpsi_channel_names.get(
            dataset_label, f"JPsi_{dataset_label}" if dataset_label else "JPsi"
        ),
        "JPsi_data": histograms["JPsi_data"],
        "JPsi_mc": histograms["JPsi_mc"],
        "JPsi_syst": histograms["nominal_muonScaleSyst_responseWeights"],
        "JPsi_resolution_syst": histograms[
            "nominal_muonResolutionSyst_responseWeights"
        ],
        "JPsi_pixel_syst": histograms["nominal_pixelMultiplicitySyst"],
    }
    if "nominal_pixelMultiplicityStat" in histograms:
        channel["JPsi_pixel_stat"] = histograms["nominal_pixelMultiplicityStat"]
    else:
        channel["JPsi_pixel_stat"] = None

    return channel


def load_jpsi_channels_from_calinput(path):
    datasets = load_calinput_datasets(path)

    # Backward compatibility for the original HDF5 layout with one data and
    # one MC dataset. The new layout uses one pair of datasets per trigger.
    if "jpsi_data" in datasets and "jpsi_mc" in datasets:
        histograms = {**datasets["jpsi_data"], **datasets["jpsi_mc"]}
        return [make_jpsi_channel("", histograms)]

    data_prefix = "jpsi_data_"
    mc_prefix = "jpsi_mc_"
    labels_data = {
        name[len(data_prefix) :] for name in datasets if name.startswith(data_prefix)
    }
    labels_mc = {
        name[len(mc_prefix) :] for name in datasets if name.startswith(mc_prefix)
    }
    if labels_data != labels_mc:
        raise RuntimeError(
            f"Mismatched JPsi trigger channels in {path}: "
            f"data-only={sorted(labels_data - labels_mc)}, "
            f"mc-only={sorted(labels_mc - labels_data)}"
        )
    if not labels_data:
        raise RuntimeError(f"No JPsi trigger-channel datasets found in {path}")

    channels = []
    for label in sorted(labels_data):
        histograms = {
            **datasets[f"{data_prefix}{label}"],
            **datasets[f"{mc_prefix}{label}"],
        }
        channels.append(make_jpsi_channel(label, histograms))
    return channels


def split_calinput_scale_variations(hsyst, hist_mc, scaling_factor, thresh):
    up_variations = []
    down_variations = []

    if "unc" not in hsyst.axes.name or "downUpVar" not in hsyst.axes.name:
        raise RuntimeError(
            "Expected WRemnants systematic histogram with 'unc' and 'downUpVar' axes"
        )

    # WRemnants uses downUpVar index 0 for down and 1 for up.
    down_index = 0
    up_index = 1
    if hsyst.axes["unc"].size % len(nuisance_categories):
        raise RuntimeError(
            f"Expected a multiple of {len(nuisance_categories)} uncertainty bins "
            f"in the calinput systematic histogram, found {hsyst.axes['unc'].size}"
        )
    num_eta_bins = hsyst.axes["unc"].size // len(nuisance_categories)
    variations = []
    groups = []

    for param_index, nuisance_category in enumerate(nuisance_categories):
        for i in range(num_eta_bins):
            unc_index = i * len(nuisance_categories) + param_index
            up = hsyst[{"unc": unc_index, "downUpVar": up_index}] * scaling_factor
            down = hsyst[{"unc": unc_index, "downUpVar": down_index}] * scaling_factor
            up_variations.append(trim_variations(clip_values(up, lim=args.clip), hist_mc, thresh))
            down_variations.append(trim_variations(clip_values(down, lim=args.clip), hist_mc, thresh))
            variations.append(f"{nuisance_category}{i}")
            groups.append(nuisance_category)

    return up_variations, down_variations, variations, groups


def split_calinput_resolution_variations(hsyst, hist_mc, scaling_factor, thresh):
    if "smearing_variation" not in hsyst.axes.name:
        raise RuntimeError(
            "Expected WRemnants resolution histogram with 'smearing_variation' axis"
        )

    res_categories = ["a", "b", "c", "d"]
    if hsyst.axes["smearing_variation"].size % len(res_categories):
        raise RuntimeError(
            f"Expected a multiple of {len(res_categories)} resolution variations, "
            f"found {hsyst.axes['smearing_variation'].size}"
        )
    num_eta_bins = hsyst.axes["smearing_variation"].size // len(res_categories)

    resolution_variations = []
    resolution_names = []
    resolution_groups = []
    for param_index, param_name in enumerate(res_categories):
        for ieta in range(num_eta_bins):
            variation_index = ieta * len(res_categories) + param_index
            variation = hsyst[{"smearing_variation": variation_index}] * scaling_factor
            resolution_variations.append(
                trim_variations(clip_values(variation, lim=args.clip), hist_mc, thresh)
            )
            resolution_names.append(f"res_{param_name}{ieta}")
            resolution_groups.append(f"res_{param_name}")

    return resolution_variations, resolution_names, resolution_groups


def split_calinput_pixel_multiplicity_variations(hsyst, hist_mc, scaling_factor, thresh, name_prefix="pixel_multiplicity_syst"):
    if "var" not in hsyst.axes.name:
        raise RuntimeError(
            f"Expected WRemnants pixel multiplicity histogram with 'var' axis, found {hsyst.axes.name}"
        )

    num_vars = hsyst.axes["var"].size
    variations = []
    names = []
    groups = []
    for ivar in range(num_vars):
        variation = hsyst[{"var": ivar}] * scaling_factor
        variations.append(
            trim_variations(clip_values(variation, lim=args.clip), hist_mc, thresh)
        )
        names.append(f"{name_prefix}_{ivar}")
        groups.append(name_prefix)

    return variations, names, groups


# set up the global systematics
nuisance_categories = ["A", "e", "M"]

print(args.systematicType)
writer = tensorwriter.TensorWriter(
    sparse=args.sparse,
    systematic_type=args.systematicType,
)

# initialize a global data covariance matrix (flattened for independent errors)
cov_data = np.array([])

# ### add D0

# D0_file = 'histsD0.pkl.lz4'
# with (lz4.frame.open(D0_file, "r")) as openfile:
#     resultdict = pickle.load(openfile)
# print(resultdict)
# hist_mc = resultdict["hD0_nom"]*0.16/2
# hist_mc_pt_up = resultdict["hD0_pt_up"]*0.16/2
# hist_mc_pt_down = resultdict["hD0_pt_down"]*0.16/2
# hist_mc_E_up = resultdict["hD0_E_up"]*0.16/2
# hist_mc_E_down = resultdict["hD0_E_down"]*0.16/2

# print(hist_mc)

# writer.add_channel(hist_mc.axes, 'D0')
# writer.add_process(hist_mc, f"sig_D0", 'D0', signal=False)
# writer.add_data(hist_mc, 'D0')

# # Axis 0 is etaK by how you loaded it (['etaK','mRK','D0mass'])
# eta_axis = hist_mc.axes[0]
# n_eta = eta_axis.size  # number of regular bins (no flow)

# nom = hist_mc.view(flow=False)
# A_up  = hist_mc_pt_up.view(flow=False)
# A_dn  = hist_mc_pt_down.view(flow=False)
# e_up  = hist_mc_E_up.view(flow=False)
# e_dn  = hist_mc_E_down.view(flow=False)

# variations_up = []
# variations_down = []

# # A variations
# for i in range(n_eta):
#     # --- UP: nominal everywhere, but eta-bin i from UP ---
#     h_up_i = hist_mc.copy()
#     v_up_i = h_up_i.view(flow=False)
#     v_up_i[i, :, :] = A_up[i, :, :]
#     variations_up.append(h_up_i)

#     # --- DOWN: nominal everywhere, but eta-bin i from DOWN ---
#     h_dn_i = hist_mc.copy()
#     v_dn_i = h_dn_i.view(flow=False)
#     v_dn_i[i, :, :] = A_dn[i, :, :]
#     variations_down.append(h_dn_i)

# # e variations
# for i in range(n_eta):
#     # --- UP: nominal everywhere, but eta-bin i from UP ---
#     h_up_i = hist_mc.copy()
#     v_up_i = h_up_i.view(flow=False)
#     v_up_i[i, :, :] = e_up[i, :, :]
#     variations_up.append(h_up_i)

#     # --- DOWN: nominal everywhere, but eta-bin i from DOWN ---
#     h_dn_i = hist_mc.copy()
#     v_dn_i = h_dn_i.view(flow=False)
#     v_dn_i[i, :, :] = e_dn[i, :, :]
#     variations_down.append(h_dn_i)

# for (up_var, down_var, var, group) in zip(variations_up, variations_down, variations, groups):
#     #print(up_var.name,var, group)
#     writer.add_systematic(
#                 [up_var, down_var],
#                 var,
#                 f"sig_D0",
#                 "D0",
#                 groups=[group],
#                 #symmetrize="none",
#                 constrained=True #should be True usually
#             )

for resultdict in load_jpsi_channels_from_calinput(args.calinputHdf5):
    sample = resultdict["name"]
    if sample not in args.jpsiChannels:
        print(f"skipping channel {sample}")
        continue
    print(f"processing {sample}")
    hist_mc = resultdict["JPsi_mc"]

    print(np.any(hist_mc.values() <= 0))
    hist_data = resultdict["JPsi_data"]
    scaling_factor = np.sum(hist_data.values()) / np.sum(hist_mc.values())
    signal_active_cells = populated_cells(hist_mc, event_thresh=args.mcEventThresh)
    print("signal active cells:", np.sum(signal_active_cells))
    print(scaling_factor, sample)
    hist_mc = hist_mc * scaling_factor
    clip_values(hist_mc, lim=args.clip, variances=True)
    zero_inactive_cells(hist_mc, signal_active_cells)
    if args.mcEventThresh < 0:
        fit_active_cells = signal_active_cells | populated_cells(hist_data)
    else:
        fit_active_cells = signal_active_cells
        zero_inactive_cells(hist_data, fit_active_cells)

    # placing a hack here to tell me which 4d bins to make rabbit selections for
    # s = ""
    # for eta1 in range(24):
    #     for eta2 in range(24):
    #         for pt1 in range(4):
    #             for pt2 in range(4):
    #                 if fit_active_cells[eta1, eta2, pt1, pt2]:
    #                     s += f"-m Select JPsi_prompt eta1:{eta1},eta2:{eta2},pt1:{pt1},pt2:{pt2} "
    # print(s)
    # exit()
    up_variations = []
    down_variations = []
    thresh = 1e-5
    up_variations, down_variations, variations, groups = (
        split_calinput_scale_variations(
            resultdict["JPsi_syst"], hist_mc, scaling_factor, thresh
        )
    )
    resolution_variations, resolution_names, resolution_groups = (
        split_calinput_resolution_variations(
            resultdict["JPsi_resolution_syst"], hist_mc, scaling_factor, thresh
        )
    )
    pixel_variations, pixel_names, pixel_groups = (
        split_calinput_pixel_multiplicity_variations(
            resultdict["JPsi_pixel_syst"], hist_mc, scaling_factor, thresh, "pixel_multiplicity_syst"
        )
    )
    if resultdict.get("JPsi_pixel_stat") is not None:
        pixel_stat_variations, pixel_stat_names, pixel_stat_groups = (
            split_calinput_pixel_multiplicity_variations(
                resultdict["JPsi_pixel_stat"], hist_mc, scaling_factor, thresh, "pixel_multiplicity_stat"
            )
        )
    else:
        pixel_stat_variations, pixel_stat_names, pixel_stat_groups = [], [], []

    writer.add_channel(hist_mc.axes, sample)
    writer.add_process(
        hist_mc, f"sig_{sample}", sample, signal=True
    )  # singal=True means cross sections float
    background_hist = make_active_background(hist_mc, fit_active_cells)
    writer.add_process(background_hist, f"background_{sample}", sample, signal=False)
    writer.add_data(hist_data, sample)
    # writer.add_norm_systematic("norm", [f"sig_{sample}"], sample, 1.02)

    # writer.add_lnN_systematic("luminosity", [f"sig_{sample}"], sample, 1.05)

    for up_var, down_var, var, group in zip(
        up_variations, down_variations, variations, groups
    ):
        print(var)
        writer.add_systematic(
            [up_var, down_var],
            var,
            f"sig_{sample}",
            sample,
            groups=[group],
            # symmetrize="none",
            constrained=True,  # should be True usually
        )

    for variation, var, group in zip(
        resolution_variations, resolution_names, resolution_groups
    ):
        print(var)
        writer.add_systematic(
            variation,
            var,
            f"sig_{sample}",
            sample,
            groups=[group],
            constrained=True,
        )

    for variation, var, group in zip(
        pixel_variations, pixel_names, pixel_groups
    ):
        print(var)
        writer.add_systematic(
            variation,
            var,
            f"sig_{sample}",
            sample,
            groups=[group],
            constrained=True,
        )

    for variation, var, group in zip(
        pixel_stat_variations, pixel_stat_names, pixel_stat_groups
    ):
        print(var)
        writer.add_systematic(
            variation,
            var,
            f"sig_{sample}",
            sample,
            groups=[group],
            constrained=True,
        )

directory = args.output
if directory == "":
    directory = "./"
filename = args.outname
if args.postfix:
    filename += f"_{args.postfix}"
writer.write(outfolder=directory, outfilename=filename)
