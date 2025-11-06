import os
import pickle
import re

import lz4.frame

from utilities import common, parsing
from wremnants import theory_corrections
from wums import boostHistHelpers as hh
from wums import logging, output_tools


def corr_name(corrf):
    if not corrf.endswith(".pkl.lz4"):
        raise ValueError(f"File {corrf} is not a lz4 compressed pickle file")

    match = re.match(r"(.*)Corr[Z|W]\.pkl\.lz4", os.path.basename(corrf))
    return match[1]


def parse_args():
    parser = parsing.base_parser()
    parser.add_argument(
        "-c",
        "--refCorr",
        type=str,
        help="Reference corr (usually a correction to another prediction+unc.)",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--rescaleCorr",
        type=str,
        help="Correction to rescale to",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--newCorrName",
        type=str,
        help="Name of the new correction",
        required=True,
    )
    parser.add_argument(
        "--outpath",
        type=str,
        default=f"{common.data_dir}/TheoryCorrections",
        help="Output path",
    )
    parser.add_argument(
        "--smooth",
        action="store_true",
        help="Apply spline-based smoothing to correction",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    ref = pickle.load(lz4.frame.open(args.refCorr))
    rescale = pickle.load(lz4.frame.open(args.rescaleCorr))

    proc = "Z" if "CorrZ" in args.refCorr else "W"

    refcorr = ref[proc][corr_name(args.refCorr) + "_minnlo_ratio"]

    ref_hist = ref[proc][corr_name(args.refCorr) + "_hist"]
    rescale_hist = rescale[proc][corr_name(args.rescaleCorr) + "_hist"]

    if ref_hist.shape[:-1] != rescale_hist.shape[:-1]:
        raise ValueError(
            f"Histograms {ref_hist.shape} and {rescale_hist.shape} have different shapes"
        )

    ratio = hh.divideHists(rescale_hist[{"vars": 0}], ref_hist[{"vars": 0}], flow=False)

    new_corr = hh.multiplyHists(refcorr, ratio, flow=False)
    if args.smooth:
        new_corr = theory_corrections.smooth_theory_corr(
            new_corr, ref[proc]["minnlo_ref_hist"], rescale_hist
        )

    output_dict = {
        args.newCorrName + "_minnlo_ratio": new_corr,
        **ref[proc],
        **rescale[proc],
    }
    meta_dict = {
        "ref_file": ref["file_meta_data"],
        "rescale_file": rescale["file_meta_data"],
    }
    outfile = f"{args.outpath}/{args.newCorrName}Corr{proc}.pkl.lz4"

    output_tools.write_lz4_pkl_output(
        outfile, proc, output_dict, common.base_dir, args, meta_dict
    )


if __name__ == "__main__":
    main()
