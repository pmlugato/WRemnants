#!/usr/bin/env bash
# The prefit/postfit readout the fit script does not produce.
#
# `btojpsik_cvhnano_fit.sh` plots the inclusive mass prefit and postfit and the
# per-cell mass POSTFIT, but not the per-cell PREFIT -- and the per-cell prefit
# is the one that says where the seed was wrong before the fit moved anything,
# which is the question when signal normalisations come back at 0.003 and 230.
# It also does not produce the per-cell uncertainty decomposition.
set -uo pipefail

WREM_BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
MZ=${MZ:-/ceph/submit/data/user/p/pmlugato/mz}
IMAGE=${IMAGE:-/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest}
CH=btojpsik_cvhnano
MASS=jointCvh_mass

run() {
  singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc \
    "cd '$WREM_BASE' && source setup.sh >/dev/null 2>&1 && $1"
}

# postfix -> cell axes present in that configuration
do_one() {
  local postfix=$1; shift
  local axes="$*"
  local pr="$MZ/fitresults/cvhnano_${postfix}_plotready.hdf5"
  local out="$MZ/plots/cvhnano_${postfix}"
  [[ -f $pr ]] || { echo "skip $postfix (no plotready)"; return; }
  mkdir -p "$out"
  echo "=== $postfix"

  # Per-cell prefit. -m must be named or the script is a silent no-op.
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists.py '$pr' -o '$out' \
       --prefit --dataName Data -m Project $CH $axes $MASS \
       --selectionAxes $axes" 2>&1 | tail -2

  # Per-cell uncertainty decomposition: which nuisance owns the band, per cell.
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists_uncertainties.py '$pr' \
       -o '$out' -m Project $CH $axes $MASS --selectionAxes $axes" 2>&1 | tail -2

  # Inclusive-mass uncertainty decomposition. No --prefit: the uncertainties
  # script has no such option (it errors out), unlike rabbit_plot_hists.py.
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists_uncertainties.py '$pr' \
       -o '$out' -m Project $CH $MASS" 2>&1 | tail -2
}

CELLS3="fit_bachelor_pt fit_bachelor_eta fit_bachelor_charge"
CELLS2="fit_bachelor_pt fit_bachelor_charge"

do_one full_exp_asimov_eta12_aem_full            $CELLS3
do_one full_exp_data_eta12_normgrouped_aem_full  $CELLS3
do_one full_exp_asimov_eta12_aem_nobbb_full      $CELLS3
do_one full_exp_asimov_eta6_aem_full             $CELLS3
do_one ptcharge_exp_data_eta12_aem_full          $CELLS2

echo done
