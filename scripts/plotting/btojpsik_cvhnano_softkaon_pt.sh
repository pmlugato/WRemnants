#!/bin/bash
# Reco kaon pT spectrum with the pT floor OPEN, through the central
# makeDataMCStackPlot.py.
#
# Input must be a histmaker output run with --bachelorMinPt 0.1 (the default is
# 1.0, which is the fit selection and would show nothing below the cut).
# The `fitbins` histogram is the one that carries fit_bachelor_pt; the standalone
# nominal_BuJpsiK_kaonPt diagnostic was not booked in the soft-kaon runs, so the
# pT axis is reached by projecting fitbins and integrating the mass window with
# --presel. Eta and charge are summed by the projection.
#
# Absolute MC normalisation is the placeholder inclusive b -> J/psi X cross
# section times lumi, so only shapes are comparable here -- the same caveat the
# rest of btojpsik_cvhnano_plots.sh carries. Raw yields are in
# scripts/btojpsik/kaon_pt_spectrum_fig.py.
#
# Usage: btojpsik_cvhnano_softkaon_pt.sh <histmaker.hdf5> [outfolder] [outpath]
#   outpath: base directory for the plots              [$HOME/public_html/mz]

set -euo pipefail

INFILE=${1:?usage: $0 <histmaker.hdf5> [outfolder] [outpath]}
OUTFOLDER=${2:-btojpsik_softkaon_pt}
OUTPATH=${3:-$HOME/public_html/mz}

WREM_BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
PLOTTER=$WREM_BASE/scripts/plotting/makeDataMCStackPlot.py

# The fit mass window, integrated away so the pT axis is the only one left.
MASS_PRESEL="jointCvh_mass=5.05,5.45"
COMMON=(-n fitbins --hists fit_bachelor_pt
        --presel "$MASS_PRESEL"
        --noRatio --legCols 1 --cmsDecor " "
        --vertLineEdges 1.0
        --customXlabel 'reco kaon $p_{T}$ [GeV]'
        -o "$OUTPATH")

# 1. The real-b stack, simulation only. NoBAncestor excluded: it is the
#    combinatorial diagnostic, not a template.
for scale in lin log; do
  extra=(--yscale 1.45); [[ $scale == log ]] && extra=(--logy --yscale 200)
  echo "--- stack, $scale"
  python "$PLOTTER" "$INFILE" "${COMMON[@]}" \
    --noData --excludeProcs NoBAncestor \
    "${extra[@]}" -f "${OUTFOLDER}_${scale}"
done

# 2. The same stack WITH the no-b-ancestor category, which is what makes the
#    soft region expensive: it is 98% of the simulated candidates below 1 GeV.
for scale in lin log; do
  extra=(--yscale 1.45); [[ $scale == log ]] && extra=(--logy --yscale 200)
  echo "--- stack with combinatorial, $scale"
  python "$PLOTTER" "$INFILE" "${COMMON[@]}" \
    --noData \
    "${extra[@]}" -f "${OUTFOLDER}_residual_${scale}"
done

# 3. Against the data. --normToData, which the mass plots in
#    btojpsik_cvhnano_plots.sh deliberately refuse, is legitimate HERE and only
#    here: the NoBAncestor group is kept, so the stack contains a combinatorial
#    component and the comparison is of the kaon pT SHAPE, not of a yield. The
#    inclusive sample still has no prompt-charmonium or fake-muon combinatorial,
#    so a residual disagreement at low pT is expected and is the point of the
#    ratio panel.
for scale in lin log; do
  extra=(--yscale 1.45); [[ $scale == log ]] && extra=(--logy --yscale 500)
  echo "--- data vs simulation shape, $scale"
  python "$PLOTTER" "$INFILE" -n fitbins --hists fit_bachelor_pt \
    --presel "$MASS_PRESEL" --legCols 1 --cmsDecor " " \
    --customXlabel 'reco kaon $p_{T}$ [GeV]' -o "$OUTPATH" \
    --vertLineEdges 1.0 --normToData --rrange 0.5 1.5 \
    --legSize 12 --extraTextLoc 0.05 0.58 \
    "${extra[@]}" -f "${OUTFOLDER}_datamc_${scale}"
done

echo "wrote to $OUTPATH/${OUTFOLDER}_{lin,log,residual_lin,residual_log,datamc_lin,datamc_log}"
