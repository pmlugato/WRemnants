#!/bin/bash
# N-1 plots for the bachelor variables, with and without the kaon pT cut.
#
# Unlike btojpsik_cvhnano_cutvars.sh -- which runs off a fully loosened
# histmaker pass and therefore shows PRE-SELECTION distributions -- both inputs
# here are at the NOMINAL working point. Gates at 0.01 / 1e6 / 200 / 1.0, muon
# pT > 3, mass [5.0, 5.5], nLegsRefit >= 3. The only thing that differs between
# the two is the kaon pT floor, so each plot is genuinely "everything applied
# except the cut on the variable being drawn".
#
#   kpt0 -- bachelorMinPt 0.   For bachelor_pt this IS the N-1 plot; for the
#           kaon-ID variables it is N-1 in pT as well, which is the configuration
#           the section 10.3 study was done in.
#   kpt1 -- bachelorMinPt 1.0. The nominal selection. For bachelor_pt this is the
#           N plot (truncated at the cut, shown for comparison); for the kaon-ID
#           variables it is the working point a threshold would actually be
#           chosen at.
#
# The kaon |eta| cut stays ON at 2.4 in both. Its axis is 48 bins over
# [-2.4, 2.4], so relaxing the cut would put the newly admitted candidates in
# the flow bins where they cannot be seen -- an N-1 plot in eta is visually
# identical to the N plot.
#
# Usage: btojpsik_cvhnano_n1_bachelor.sh <kpt0.hdf5> <kpt1.hdf5> [outfolder] [outpath]
#   outpath: base directory for the plots              [$HOME/public_html/mz]

set -euo pipefail

F0=${1:?usage: $0 <kpt0.hdf5> <kpt1.hdf5> [outfolder] [outpath]}
F1=${2:?usage: $0 <kpt0.hdf5> <kpt1.hdf5> [outfolder] [outpath]}
OUTFOLDER=${3:-cutvars}
OUTPATH=${4:-$HOME/public_html/mz}

WREM_BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
PLOTTER=$WREM_BASE/scripts/plotting/makeDataMCStackPlot.py

# hist:axis:vlines:extra
VARS=(
  # The B mass, which is what the cuts are ultimately protecting. Lines mark
  # the fit window [5.05, 5.45], narrower than the [5.0, 5.5] selection.
  "nominal_BuJpsiK_jointCvhMass:jointCvh_mass:5.05 5.45:--yscale 1.4"
  "nominal_BuJpsiK_jointCvhMass:jointCvh_mass:5.05 5.45:--logy --yscale 300"
  "nominal_BuJpsiK_kaonPt:bachelor_pt:1.0:--logy --yscale 300"
  "nominal_BuJpsiK_kaonEta:bachelor_eta:-:--yscale 1.4"
  "nominal_cand_bachelorNValidHits:bachelor_nValidHits:-:--logy --yscale 300"
  "nominal_cand_bachelorD0:bachelor_d0:-:--logy --yscale 300"
  # The kaon ID as actually measured: |d0|/sigma > 1.5 took the exponential
  # background chi2/ndf from 41.2 to 25.4 at 86% signal efficiency
  # (add-cvhnano-bachelor-quality-cuts 10.3) and is in no fit on disk.
  "nominal_cand_bachelorD0Sig:bachelor_d0Sig:-1.5 1.5:--logy --yscale 300"
)

run_set() {
  local infile=$1 sub=$2
  for entry in "${VARS[@]}"; do
    local hname axis vl extra vlarg
    IFS=':' read -r hname axis vl extra <<< "$entry"
    vlarg=()
    [[ "$vl" != "-" ]] && vlarg=(--vertLineEdges $vl)
    echo "--- $sub  $axis"
    # shellcheck disable=SC2086
    python "$PLOTTER" "$infile" \
      -n "$hname" --hists "$axis" \
      --normToData --ratioToData --rrange 0.5 1.5 \
      --legCols 1 --legSize 15 --cmsDecor " " \
      --noLowerLeg --stackOrder NoBAncestor --extraTextLoc 0.04 0.28 \
      "${vlarg[@]}" $extra \
      -o "$OUTPATH" -f "$OUTFOLDER/$sub" || echo "    FAILED: $axis"
  done
}

run_set "$F0" n1_nokptcut
run_set "$F1" n1_kpt1

echo "wrote to $OUTPATH/$OUTFOLDER/{n1_nokptcut,n1_kpt1}"
