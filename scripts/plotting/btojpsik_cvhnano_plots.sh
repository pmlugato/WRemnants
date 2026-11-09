#!/bin/bash
# Prefit plots for the CVH NanoAOD B+ -> J/psi K histmaker output.
#
# Two stages, and the split is the point.
#
# Stage 1 -- SIMULATION ONLY, via the central makeDataMCStackPlot.py. The
# genCategory axis is sliced into stacked components by
# wremnants/postprocessing/datagroups/datagroups_cvhnano.py, so the standard
# tool needs nothing beyond the labels and colours in
# wremnants/utilities/styles/styles.py.
#
#   --noData, and NOT --normToData, for a physics reason. Every simulated component is a real b-hadron decay; the
#   simulation contains none of the prompt-charmonium or fake-muon combinatorial
#   that data also has, and the combinatorial it does contain is excluded from
#   the groups on purpose. The stack is therefore ~12% of the data in the peak
#   window and ~0.5% outside it. Normalising it to the data yield would scale the
#   signal up by roughly eight to cover background the sample does not have, and
#   the ratio panel would then read as agreement by construction.
#
# Stage 2 -- DATA, via btojpsik_cvhnano_prefit.py, which fits the signal
# template plus an analytic background to the data so that the signal
# normalisation is measured rather than imposed. Mass variables only.
#
# --legCols 1 rather than the default 2: a two-column legend puts its left
# column on top of the CMS label. --noRatio because with --noData there is
# nothing to take a ratio to, and the tool still draws the empty panel.
#
# Usage: btojpsik_cvhnano_plots.sh <histmaker.hdf5> [outfolder] [outpath]
#   outpath: base directory for the plots              [$HOME/public_html/mz]

set -euo pipefail

INFILE=${1:?usage: $0 <histmaker.hdf5> [outfolder] [outpath]}
OUTFOLDER=${2:-btojpsik_cvhnano}
OUTPATH=${3:-$HOME/public_html/mz}

WREM_BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
PLOTTER=$WREM_BASE/scripts/plotting/makeDataMCStackPlot.py
PREFIT=$WREM_BASE/scripts/plotting/btojpsik_cvhnano_prefit.py

# hist-name:axis-name pairs. The histmaker writes one histogram per variable,
# each carrying its own axis plus genCategory.
VARS=(
  # B candidate, one entry per fit arm
  "nominal_BuJpsiK_jointCvhMass:jointCvh_mass"
  "nominal_BuJpsiK_kvfCvhMass:kvfCvh_mass"
  "nominal_BuJpsiK_kvfRawMass:kvfRaw_mass"
  "nominal_BuJpsiK_mass:cand_mass"
  "nominal_BuJpsiK_jointCvhMassErr:jointCvh_massErr"
  "nominal_BuJpsiK_jointCvhPt:jointCvh_pt"
  # J/psi. jointCvhJpsiMass is the constrained parameter and is deliberately
  # included to show that it carries no information; the raw-track rebuild
  # and the KVF fits below are the unconstrained dimuon masses.
  "nominal_cand_dimuonRawMass:dimuon_raw_mass"
  "nominal_cand_dimuonRawPt:dimuon_raw_pt"
  "nominal_BuJpsiK_jointCvhJpsiMass:jointCvh_jpsi_mass"
  # Unconstrained dimuon mass from a real vertex fit. This is
  # the muon momentum-scale handle; the two above are a raw-track rebuild and a
  # constrained parameter respectively.
  "nominal_BuJpsiK_kvfCvhDimuonMass:kvfCvh_dimuon_mass"
  "nominal_BuJpsiK_kvfRawDimuonMass:kvfRaw_dimuon_mass"
  # Legs
  "nominal_cand_mu0Pt:mu0_pt"
  "nominal_cand_mu1Pt:mu1_pt"
  "nominal_cand_mu0Eta:mu0_eta"
  "nominal_BuJpsiK_kaonPt:bachelor_pt"
  "nominal_BuJpsiK_kaonEta:bachelor_eta"
  # Impact parameters, referenced to the primary vertex.
  "nominal_cand_bachelorD0:bachelor_d0"
  "nominal_cand_bachelorDzPV:bachelor_dzPV"
  "nominal_cand_bachelorNValidHits:bachelor_nValidHits"
  # Impact-parameter significance and track quality.
  "nominal_cand_bachelorD0Sig:bachelor_d0Sig"
  "nominal_cand_bachelorRelPtErr:bachelor_relPtErr"
  "nominal_cand_bachelorNValidPixelHits:bachelor_nValidPixelHits"
  "nominal_cand_bachelorTrackerLayers:bachelor_trackerLayers"
  # Vertex and displacement
  "nominal_BuJpsiK_jointCvhSl3d:jointCvh_sl3d"
  "nominal_BuJpsiK_jointCvhLxy:jointCvh_lxy"
  "nominal_BuJpsiK_jointCvhAlphaBS:jointCvh_alphaBS"
  # 3D pointing angle to the PV.
  "nominal_BuJpsiK_jointCvhAlpha3dPV:jointCvh_alpha3dPV"
  "nominal_BuJpsiK_kvfCvhVtxProb:kvfCvh_vtxProb"
  # Fit diagnostics and event quantities
  "nominal_BuJpsiK_jointCvhChisq:jointCvh_chisq"
  "nominal_BuJpsiK_jointCvhNiter:jointCvh_nIter"
  "nominal_BuJpsiK_nLegsRefit:nLegsRefit"
  "nominal_nBuJpsiK:nCand"
  "nominal_nPV:nPV"
  "nominal_nTrack:nTrack"
  # True pileup, the variable the reweighting is keyed on; nPV above is then
  # a test of it.
  "nominal_Pileup_nTrueInt:nTrueInt"
)

for entry in "${VARS[@]}"; do
  hname=${entry%%:*}
  axis=${entry##*:}
  echo "--- $axis"
  python "$PLOTTER" "$INFILE" \
    -n "$hname" --hists "$axis" \
    --noData --noRatio --yscale 1.45 --legCols 1 \
    --excludeProcs NoBAncestor \
    --cmsDecor " " \
    -o "$OUTPATH" -f "$OUTFOLDER"
done

# The one plot that DOES include the combinatorial diagnostic group, because its
# whole purpose is to show what is left over once the real-b decays are named.
# Log scale: the residual is ~78% of the peak window, so on a linear axis the
# real-b components are invisible next to it.
echo "--- residual: the same stack WITH the no-b-ancestor category"
for entry in "nominal_BuJpsiK_jointCvhMass:jointCvh_mass" \
             "nominal_BuJpsiK_kvfCvhMass:kvfCvh_mass"; do
  hname=${entry%%:*}
  axis=${entry##*:}
  python "$PLOTTER" "$INFILE" \
    -n "$hname" --hists "$axis" \
    --noData --noRatio --yscale 1.45 --legCols 1 \
    --cmsDecor " " --logy \
    -o "$OUTPATH" -f "${OUTFOLDER}_residual"
done

echo "--- data-side prefit (mass variables)"
python "$PREFIT" "$INFILE" -o "$OUTPATH" --outfolder "${OUTFOLDER}_prefit"

echo "wrote to $OUTPATH/$OUTFOLDER and $OUTPATH/${OUTFOLDER}_prefit"
