#!/bin/bash
# Every variable the CVH NanoAOD B+ histmaker cuts on, drawn with the cut on it,
# through the central makeDataMCStackPlot.py.
#
# Input must be a LOOSENED histmaker run -- gates at 1e9, --muonMinPt 0
# --bachelorMinPt 0 --minLegsRefit 0, --allaxes. A variable plotted after its
# own cut is truncated exactly where the interesting part is. The arm-ok flag
# and the sentinel rejection stay on even in the loose run: a failed fit writes
# -99, not NaN, so admitting sentinels corrupts every distribution on the page.
#
# Two differences from btojpsik_cvhnano_plots.sh, both deliberate.
#
#   NoBAncestor is INCLUDED here. That script excludes it because a template
#   built from one third of the combinatorial would look like a measurement of
#   all of it. But the question these plots ask is whether a cut separates
#   signal from combinatorial, and the simulation's combinatorial is the only
#   handle on that, so hiding it would remove the subject of the plot.
#
#   --normToData, which that script argues against for the mass stack and for
#   the same reason: the simulation carries no prompt-charmonium or fake-muon
#   combinatorial, so it is ~12% of data in the peak window and normalising
#   inflates the signal to cover background the sample does not have. Here the
#   comparison is explicitly of SHAPE -- does the simulation describe the
#   variable we are cutting on -- so the ratio panel has to read against 1. The
#   absolute deficit is a separate plot, not this one.
#
# --vertLineEdges draws the threshold on both panels. It needs the ratio panel,
# so --noRatio and a vertical line cannot be combined.
#
# Usage: btojpsik_cvhnano_cutvars.sh <loose_histmaker.hdf5> [outfolder] [group]
#   group: all (default) | active | off | jpsi

set -euo pipefail

INFILE=${1:?usage: $0 <loose_histmaker.hdf5> [outfolder]}
OUTFOLDER=${2:-cutvars_260908}
GROUP=${3:-all}
OUTPATH=${OUTPATH:-$HOME/public_html/mz}

WREM_BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
PLOTTER=$WREM_BASE/scripts/plotting/makeDataMCStackPlot.py

# hist:axis:vlines:extra
#   vlines "-"  -> no line (threshold off the plotted range, or a bare flag)
#   extra       -> per-variable flags. A one-column legend needs a lot of
#                  headroom on a log axis, hence the large --yscale there.
ACTIVE=(
  # edmRef bin edges are log-spaced over 1e-11..1e2, so on a linear x axis the
  # whole distribution piles up at zero and the gate at 0.01 is invisible.
  "nominal_BuJpsiK_jointCvhEdmRef:jointCvh_edmRef:0.01:--logy --logx --yscale 300"
  "nominal_BuJpsiK_jointCvhChisq:jointCvh_chisq:-:--logy --yscale 300"
  "nominal_BuJpsiK_jointCvhNdof:jointCvh_ndof:-:--logy --yscale 300"
  "nominal_BuJpsiK_jointCvhMassErr:jointCvh_massErr:-:--logy --yscale 300"
  "nominal_cand_mu0Pt:mu0_pt:3.0:--yscale 1.4"
  "nominal_cand_mu1Pt:mu1_pt:3.0:--yscale 1.4"
  "nominal_cand_mu0Eta:mu0_eta:-:--yscale 1.4"
  "nominal_cand_mu1Eta:mu1_eta:-:--yscale 1.4"
  "nominal_BuJpsiK_kaonPt:bachelor_pt:1.0:--logy --yscale 300"
  "nominal_BuJpsiK_kaonEta:bachelor_eta:-:--yscale 1.4"
  "nominal_BuJpsiK_jointCvhMass:jointCvh_mass:5.05 5.45:--yscale 1.4"
  "nominal_BuJpsiK_nLegsRefit:nLegsRefit:-:--logy --yscale 300"
)

# Declared, at their no-op defaults. Thresholds are the values the studies used:
# |d0|/sigma > 1.5 is the kaon ID from add-cvhnano-bachelor-quality-cuts 10.3;
# the dimuon and B-vertex values are 2018's, from section 10.2.
OFF=(
  "nominal_cand_bachelorD0Sig:bachelor_d0Sig:-1.5 1.5:--logy --yscale 300"
  "nominal_cand_bachelorD0:bachelor_d0:-:--logy --yscale 300"
  "nominal_cand_bachelorDzPVSig:bachelor_dzPVSig:-:--logy --yscale 300"
  "nominal_cand_bachelorDedxHarmonic2:bachelor_dedx:-:--logy --yscale 300"
  "nominal_cand_bachelorNValidHits:bachelor_nValidHits:-:--logy --yscale 300"
  "nominal_cand_bachelorNormChi2:bachelor_normChi2:-:--yscale 1.4"
  "nominal_cand_bachelorRelPtErr:bachelor_relPtErr:-:--yscale 1.4"
  "nominal_BuJpsiK_jointCvhNiter:jointCvh_nIter:-:--logy --yscale 300"
  "nominal_BuJpsiK_jointCvhSl3d:jointCvh_sl3d:-:--logy --yscale 300"
  "nominal_BuJpsiK_jointCvhAlphaBS:jointCvh_alphaBS:-:--logy --yscale 300"
  "nominal_BuJpsiK_jointCvhAlpha3dPV:jointCvh_alpha3dPV:-:--logy --yscale 300"
  "nominal_BuJpsiK_kvfCvhVtxProb:kvfCvh_vtxProb:0.1 0.3:--logy --yscale 300"
  "nominal_BuJpsiK_kvfCvhDimuonVtxProb:kvfCvh_dimuon_vtxProb:0.1:--logy --yscale 300"
  "nominal_BuJpsiK_kvfCvhDimuonSl3d:kvfCvh_dimuon_sl3d:-:--logy --yscale 300"
  "nominal_BuJpsiK_kvfCvhDimuonAlphaBS:kvfCvh_dimuon_alphaBS:0.4:--logy --yscale 300"
  "nominal_BuJpsiK_jpsiPt:jpsi_raw_pt:7.0:--yscale 1.4"
)

# Not cut on. The J/psi kinematics as the joint N-body CVH fit returns them
# (BuJpsiK_jointCvhJpsiPt/Eta), i.e. with the J/psi mass constraint and the
# B vertex applied, as opposed to jpsi_raw_pt above, which is the pre-fit
# dimuon. Same 7 GeV marker as jpsi_raw_pt, 2018's J/psi pT floor, for
# comparison only: nothing in the histmaker cuts on it.
JPSI=(
  "nominal_BuJpsiK_jointCvhJpsiPt:jointCvh_jpsi_pt:7.0:--yscale 1.4"
  "nominal_BuJpsiK_jointCvhJpsiPt:jointCvh_jpsi_pt:7.0:--logy --yscale 300"
  "nominal_BuJpsiK_jointCvhJpsiEta:jointCvh_jpsi_eta:-:--yscale 1.4"
)

run_one() {
  local entry=$1 sub=$2
  local hname axis vl extra
  IFS=':' read -r hname axis vl extra <<< "$entry"
  local vlarg=()
  [[ "$vl" != "-" ]] && vlarg=(--vertLineEdges $vl)
  echo "--- $sub  $axis"
  # shellcheck disable=SC2086
  python "$PLOTTER" "$INFILE" \
    -n "$hname" --hists "$axis" \
    --normToData --ratioToData --rrange 0.5 1.5 \
    --legCols 1 --legSize 15 --cmsDecor " " \
    --noLowerLeg --stackOrder NoBAncestor --extraTextLoc 0.04 0.28 \
    "${vlarg[@]}" $extra \
    -o "$OUTPATH" -f "$OUTFOLDER/$sub" || echo "    FAILED: $axis"
}

case $GROUP in
  all|active) for e in "${ACTIVE[@]}"; do run_one "$e" active; done ;;&
  all|off)    for e in "${OFF[@]}";    do run_one "$e" declared_off; done ;;&
  all|jpsi)   for e in "${JPSI[@]}";   do run_one "$e" jpsi_jointCvh; done ;;
  *) echo "unknown group '$GROUP' (all|active|off|jpsi)" >&2; exit 2 ;;
esac

echo "wrote to $OUTPATH/$OUTFOLDER"
