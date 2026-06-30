import argparse
import os
import sys

import h5py


def print_flush(*args, **kwargs):
    print(*args, **kwargs)
    sys.stdout.flush()


wremnants_base = os.environ.get("WREM_BASE", None)
if wremnants_base is None:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    wremnants_base = os.path.abspath(os.path.join(script_dir, "../.."))
    if os.path.exists(os.path.join(wremnants_base, "wums")):
        sys.path.insert(0, wremnants_base)

from rabbit import tensorwriter
from wremnants.utilities.io_tools import base_io

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Input HDF5 file with histograms")
parser.add_argument("-o", "--output", default="./", help="output directory")
parser.add_argument("--outname", default="my_tensor", help="output file name")
parser.add_argument(
    "--histName", default="ptll", help="Histogram name to use (default: ptll)"
)
parser.add_argument(
    "--procFilters",
    nargs="*",
    default=["Zmumu"],
    help="Processes to include (default: Zmumu)",
)
parser.add_argument(
    "--sparse", default=False, action="store_true", help="Make sparse tensor"
)
parser.add_argument(
    "--systematicType",
    choices=["log_normal", "normal"],
    default="log_normal",
    help="probability density for systematic variations",
)
args = parser.parse_args()

# Load data from HDF5 file
infile_path = os.path.abspath(args.infile)
if not os.path.exists(infile_path):
    raise FileNotFoundError(f"Input file not found: {infile_path}")

h5file = h5py.File(infile_path, "r")
results = base_io.load_results_h5py(h5file)

hist_name = args.histName
all_procs = [p for p in results.keys() if p != "meta_info"]
print_flush(f"All processes found in file: {all_procs}")

if args.procFilters:
    procs_to_use = [p for p in all_procs if any(filt in p for filt in args.procFilters)]
else:
    procs_to_use = all_procs

print_flush(f"Processes after filtering: {procs_to_use}")

# All processes are treated as MC (no data loading)
mc_procs = procs_to_use
print_flush(f"MC processes found: {mc_procs}")

# Load MC histograms
h_mc_dict = {}
for proc in mc_procs:
    if "output" in results[proc] and hist_name in results[proc]["output"]:
        h_proxy = results[proc]["output"][hist_name]
        h_mc_dict[proc] = h_proxy.get() if hasattr(h_proxy, "get") else h_proxy

# Identify signal processes before closing file
signal_procs = [p for p in mc_procs if "Zmumu" in p or "Ztautau" in p]


def _theory_corr_hist_name(results, proc, tag):
    """Match histmaker output: <base>_<generator>_Corr (e.g. ptll_..._pdfas_Corr)."""
    for name in results[proc]["output"]:
        if name.endswith("_Corr") and tag in name:
            return name
    return None


h_theory_corr_pdfas = None
h_theory_corr_pdfvars = None
h_theory_corr_pdfas_dict = {}
h_theory_corr_pdfvars_dict = {}

if signal_procs:
    first_sig_proc = signal_procs[0]
    if first_sig_proc in results and "output" in results[first_sig_proc]:
        theory_corr_hist_name_pdfas = _theory_corr_hist_name(
            results, first_sig_proc, "pdfas"
        )
        theory_corr_hist_name_pdfvars = _theory_corr_hist_name(
            results, first_sig_proc, "pdfvars"
        )
        if not theory_corr_hist_name_pdfas:
            print_flush(
                f"Warning: no pdfas theory hist in {first_sig_proc}; "
                f"keys: {list(results[first_sig_proc]['output'].keys())}"
            )
        if not theory_corr_hist_name_pdfvars:
            print_flush(f"Warning: no pdfvars theory hist in {first_sig_proc}")
        # Load pdfas correction (for alpha_s variations)
        if (
            theory_corr_hist_name_pdfas
            and theory_corr_hist_name_pdfas in results[first_sig_proc]["output"]
        ):
            h_proxy = results[first_sig_proc]["output"][theory_corr_hist_name_pdfas]
            h_theory_corr_pdfas = h_proxy.get() if hasattr(h_proxy, "get") else h_proxy
            print_flush(
                f"Found theory correction histogram for alpha_s variations: {theory_corr_hist_name_pdfas}"
            )
            print_flush(
                f"Histogram axes: {[ax.name for ax in h_theory_corr_pdfas.axes]}"
            )

            # Load for all signal processes
            for proc in signal_procs:
                if proc in results and "output" in results[proc]:
                    if theory_corr_hist_name_pdfas in results[proc]["output"]:
                        h_proxy = results[proc]["output"][theory_corr_hist_name_pdfas]
                        h_theory_corr_pdfas_dict[proc] = (
                            h_proxy.get() if hasattr(h_proxy, "get") else h_proxy
                        )

        # Load pdfvars correction (for PDF variations)
        if (
            theory_corr_hist_name_pdfvars
            and theory_corr_hist_name_pdfvars in results[first_sig_proc]["output"]
        ):
            h_proxy = results[first_sig_proc]["output"][theory_corr_hist_name_pdfvars]
            h_theory_corr_pdfvars = (
                h_proxy.get() if hasattr(h_proxy, "get") else h_proxy
            )
            print_flush(
                f"Found theory correction histogram for PDF variations: {theory_corr_hist_name_pdfvars}"
            )
            print_flush(
                f"Histogram axes: {[ax.name for ax in h_theory_corr_pdfvars.axes]}"
            )

            # Load for all signal processes
            for proc in signal_procs:
                if proc in results and "output" in results[proc]:
                    if theory_corr_hist_name_pdfvars in results[proc]["output"]:
                        h_proxy = results[proc]["output"][theory_corr_hist_name_pdfvars]
                        h_theory_corr_pdfvars_dict[proc] = (
                            h_proxy.get() if hasattr(h_proxy, "get") else h_proxy
                        )

h5file.close()

# Use first MC process as expected data (Asimov data)
if h_mc_dict:
    first_mc_proc = list(h_mc_dict.keys())[0]
    h_data = h_mc_dict[first_mc_proc].copy()
    print_flush(f"Using MC process '{first_mc_proc}' as expected data (Asimov)")
else:
    raise RuntimeError("No MC processes found to use as data")

# Handle vars axis if present
axis_names = [ax.name for ax in h_data.axes] if hasattr(h_data, "axes") else []
has_vars = "vars" in axis_names

if has_vars:
    h_data_base = h_data[{"vars": 0}]
    h_mc_base = {proc: h[{"vars": 0}] for proc, h in h_mc_dict.items()}
else:
    h_data_base = h_data
    h_mc_base = h_mc_dict

# Build tensor
writer = tensorwriter.TensorWriter(
    sparse=args.sparse,
    systematic_type=args.systematicType,
)

channel_name = "ch0"
writer.add_channel(h_data_base.axes, channel_name)
writer.add_data(h_data_base, channel_name)

background_procs = [p for p in mc_procs if p not in signal_procs]

# Add processes
for proc in signal_procs:
    if proc in h_mc_base:
        writer.add_process(h_mc_base[proc], proc, channel_name, signal=False)

for proc in background_procs:
    if proc in h_mc_base:
        writer.add_process(h_mc_base[proc], proc, channel_name, signal=False)

# Extract PDF variations from pdfvars correction histogram (noi=False, default)
if h_theory_corr_pdfvars is not None and signal_procs:
    proc_name = signal_procs[0]
    h_sig_corr_pdfvars = h_theory_corr_pdfvars_dict.get(
        proc_name, h_theory_corr_pdfvars
    )

    if "vars" in h_sig_corr_pdfvars.axes.name:
        vars_axis = h_sig_corr_pdfvars.axes["vars"]
        print_flush(
            f"Found pdfvars vars axis with entries: {[str(v) for v in vars_axis]}"
        )

        # Extract all PDF variations (exclude alpha_s variations which have "_as_" in the name)
        pdf_variations = []
        for var_name in vars_axis:
            var_str = str(var_name)
            # Skip central (index 0) and alpha_s variations
            if var_str != "central" and "_as_" not in var_str:
                pdf_variations.append(var_str)

        if pdf_variations:
            print_flush(
                f"Found {len(pdf_variations)} PDF variations: {pdf_variations[:5]}..."
            )  # Print first 5

            for var_name_up, var_name_down in zip(
                pdf_variations[1::2], pdf_variations[2::2]
            ):
                h_var_up = h_sig_corr_pdfvars[{"vars": var_name_up}]
                h_var_down = h_sig_corr_pdfvars[{"vars": var_name_down}]

                h_pdf_var_up = h_var_up.project(*h_mc_base[proc_name].axes.name)
                h_pdf_var_down = h_var_down.project(*h_mc_base[proc_name].axes.name)

                syst_name = f"{var_name_up}_{var_name_down}"

                writer.add_systematic(
                    [h_pdf_var_up, h_pdf_var_down],  # List of [up, down] histograms
                    syst_name,
                    proc_name,
                    "ch0",
                    constrained=True,
                    # mirror=False
                )

            num_pairs = len(pdf_variations) // 2
            print_flush(
                f"Added {num_pairs} PDF systematic pairs for process {proc_name} (noi=False)"
            )
        else:
            print_flush(f"Warning: No PDF variations found in pdfvars histogram")
    else:
        print_flush(
            f"Warning: pdfvars correction histogram does not have vars axis. Axes: {[ax.name for ax in h_sig_corr_pdfvars.axes]}"
        )

# Extract alphas variations from pdfas correction histogram (noi=True)
if h_theory_corr_pdfas is not None and signal_procs:
    proc_name = signal_procs[0]
    h_sig_corr_pdfas = h_theory_corr_pdfas_dict.get(proc_name, h_theory_corr_pdfas)

    if "vars" in h_sig_corr_pdfas.axes.name:
        vars_axis = h_sig_corr_pdfas.axes["vars"]
        print_flush(
            f"Found pdfas vars axis with entries: {[str(v) for v in vars_axis]}"
        )

        var_name_0120 = "pdfCT18ZNNLO_as_0120"
        var_name_0116 = "pdfCT18ZNNLO_as_0116"

        h_alphas_0120 = None
        h_alphas_0116 = None

        if var_name_0120 in vars_axis:
            h_var = h_sig_corr_pdfas[{"vars": var_name_0120}]
            h_alphas_0120 = h_var.project(*h_mc_base[proc_name].axes.name)
            print_flush(f"Found alphas 0120 variation: {var_name_0120}")

        if var_name_0116 in vars_axis:
            h_var = h_sig_corr_pdfas[{"vars": var_name_0116}]
            h_alphas_0116 = h_var.project(*h_mc_base[proc_name].axes.name)
            print_flush(f"Found alphas 0116 variation: {var_name_0116}")

        if h_alphas_0120 is not None and h_alphas_0116 is not None:
            writer.add_systematic(
                [h_alphas_0120, h_alphas_0116],
                "pdfAlphaS",
                proc_name,
                "ch0",
                constrained=False,
                noi=True,  # Only alpha_s variations have noi=True
            )
            print_flush(
                f"Added pdfAlphaS systematic for process {proc_name} (noi=True)"
            )
        else:
            print_flush(f"Warning: Could not find both alpha_s variations")
    else:
        print_flush(
            f"Warning: pdfas correction histogram does not have vars axis. Axes: {[ax.name for ax in h_sig_corr_pdfas.axes]}"
        )

# Write output
writer.write(outfolder=args.output, outfilename=args.outname)
print_flush(f"Tensor written to: {args.output}/{args.outname}.hdf5")
