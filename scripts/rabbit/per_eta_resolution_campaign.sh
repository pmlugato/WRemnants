#!/usr/bin/env bash
# Run matrix of openspec/changes/add-cvhnano-per-eta-resolution, as driver
# invocations. Every arm shares the primary binning and the pinned signal scale;
# the arms differ only in the lineshape ingredients, so each one answers one
# question of the proposal.
#
# Usage: per_eta_resolution_campaign.sh [--arms LIST] [--observed] [--hist kpi|muon_cor]
#   --arms LIST   comma-separated subset of
#                 primary,constrained,etapt,global,nores,notmpl,split,ref003,bare,ref006
#                 [all of them]. `bare` is the 2026-09-08 model (no resolution,
#                 no template) at this binning and pinned scale: the reference
#                 that says whether a data-side failure is new.
#   --observed    fit the observation instead of the Asimov dataset
#   --hist H      kpi (default) or muon_cor
#   --dryRun      pass --dryRun to the driver
#
# The pinned scale: 1,327,932 (data-fitted N_sig, aem 11.3, [5.05, 5.45]) over
# 72,790 (unscaled simulated signal in the same window on both kpiA histograms,
# measured 2026-09-13 from the writer metadata) = 18.2433.
set -uo pipefail
WREM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
HISTS=/ceph/submit/data/user/p/pmlugato/mz/histograms
ARMS="primary,constrained,etapt,global,nores,notmpl,split,ref003,bare"
OBS=""; HIST=kpi; DRY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --arms) ARMS=$2; shift 2 ;;
    # A calibration channel has nothing to blind, and rabbit's default blinding
    # multiplies every normalisation by exp(N(0, 5)) on an observed fit -- which
    # is what made every B+ data fit before 2026-09-14 stall (prefit loss 986k
    # blinded, 46k unblinded).
    --observed) OBS="--observed"; shift ;;
    --hist) HIST=$2; shift 2 ;;
    --dryRun) DRY="--dryRun"; shift ;;
    *) echo "unknown option $1" >&2; exit 2 ;;
  esac
done
SIGNAL_SCALE=18.2433
REF_YIELD=1327932
COMMON=(full --bkgModel exppoly --etaRebin 2 --ptBins 3 --normAxes fit_bachelor_eta
        --scaleVar --signalScale $SIGNAL_SCALE --refYield $REF_YIELD
        --massLow 5.05 --massHigh 5.45 --massBins 20 $OBS $DRY)
# The driver unblinds observed fits itself; the suffix keeps the provenance
# visible next to the blinded first-pass outputs of 2026-09-14.
if [[ -n $OBS ]]; then COMMON+=(--postfixSuffix "_${HIST}_unblind"); else COMMON+=(--postfixSuffix "_$HIST"); fi
H="$HISTS/btojpsik_cvhnano_kpiA_${HIST}.hdf5"
run_arm() { [[ ",$ARMS," == *",$1,"* ]] || return 0; shift
  echo; echo "########## $(date '+%F %T')  arm: $* "; echo
  "$WREM/scripts/rabbit/btojpsik_cvhnano_fit.sh" "$H" "${COMMON[@]}" "$@"; }
run_arm primary     --mcBkg peak  --resVar 0.005 --resAxes fit_bachelor_eta
run_arm constrained --mcBkg peak  --resVar 0.005 --resAxes fit_bachelor_eta --resConstrained
run_arm etapt       --mcBkg peak  --resVar 0.005 --resAxes fit_bachelor_eta,fit_bachelor_pt
run_arm global      --mcBkg peak  --resVar 0.005
run_arm nores       --mcBkg peak
run_arm notmpl      --mcBkg none  --resVar 0.005 --resAxes fit_bachelor_eta
run_arm split       --mcBkg split --resVar 0.005 --resAxes fit_bachelor_eta
run_arm ref003      --mcBkg peak  --resVar 0.003 --resAxes fit_bachelor_eta
run_arm bare        --mcBkg none
run_arm ref006      --mcBkg peak  --resVar 0.006 --resAxes fit_bachelor_eta
echo; echo "########## $(date '+%F %T')  campaign done"
