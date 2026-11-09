#!/usr/bin/env bash
# Prefit, postfit, correlation and uncertainty readout for the Asimov fit under
# BOTH background models, so the comparison isolates one thing: does giving the
# background a curvature parameter cost the A/e/M uncertainties anything?
#
# Asimov only, deliberately. The data arm does not converge (edmval 72668 on exp,
# indefinite Hessian on expquad), so a data correlation matrix is not a
# covariance. Asimov converges to edmval ~1e-28 under both models, which is
# exactly what makes this comparison meaningful: the uncertainties come from the
# Hessian at a genuine minimum, and the only difference between the two runs is
# the one parameter.
#
# Usage: cvhnano_expquad_asimov_plots.sh
set -uo pipefail

WREM_BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
MZ=${MZ:-/ceph/submit/data/user/p/pmlugato/mz}
OUT=${OUT:-$HOME/public_html/mz/cvhnano-aem/expquad-asimov}
IMAGE=${IMAGE:-/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest}
CH=btojpsik_cvhnano
MASS=jointCvh_mass

run() {
  singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc \
    "cd '$WREM_BASE' && source setup.sh >/dev/null 2>&1 && $1"
}

for MODEL in exp expquad; do
  TAG="full_${MODEL}_asimov_eta12_normgrouped_aem_full"
  FIT="$MZ/fitresults/cvhnano_${TAG}_fit.hdf5"
  PLOTREADY="$MZ/fitresults/cvhnano_${TAG}_plotready.hdf5"
  D="$OUT/$MODEL"
  if [[ ! -f $FIT ]]; then echo "missing $FIT" >&2; continue; fi
  mkdir -p "$D"
  echo "=== $MODEL : $(basename "$FIT")"

  # Prefit and postfit on the inclusive mass projection. rabbit_plot_hists.py
  # takes --prefit; rabbit_plot_hists_uncertainties.py does NOT, which is the
  # inconsistency to remember rather than rediscover.
  for WHEN in "--prefit prefit" ": postfit"; do
    set -- $WHEN; FLAG=$1; NAME=$2
    [[ $FLAG == ":" ]] && FLAG=""
    run "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists.py '$PLOTREADY' \
         -o '$D' -p '${NAME}_incl' $FLAG \
         -m Project $CH $MASS --dataName Asimov" 2>&1 | tail -2
  done

  # Uncertainty decomposition: the point of the exercise.
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists_uncertainties.py \
       '$PLOTREADY' -o '$D' -p 'unc_incl' \
       -m Project $CH $MASS" 2>&1 | tail -2

  # Correlations, sliced the same three ways as the data readout so the two are
  # comparable: A/e/M inside a group, one parameter across groups, A+e global.
  cov() { local tag=$1 lab=$2; shift 2
    run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -v 4 -o '$D' \
         --correlation --params $* --showNumbers -p '$tag' \
         --xlabel '$lab' --ylabel '$lab'" 2>&1 | tail -2
  }
  for g in 0 1 2 3; do
    cov "eta${g}_aem" "A/e/M eta$g" "scaleA${g}" "scalee${g}" "scaleM${g}"
  done
  # NB the e nuisance is "scalee", lower case, unlike scaleA and scaleM.
  for p in A e M; do
    cov "${p}_across_eta" "$p by eta" "scale${p}0" "scale${p}1" "scale${p}2" "scale${p}3"
  done
  # 8x8 without numbers -- with them it is an unreadable smear, and the pattern
  # (block diagonal or not) is the whole question.
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -v 4 -o '$D' \
       --correlation -p 'Ae_all' --params \
       scaleA0 scaleA1 scaleA2 scaleA3 scalee0 scalee1 scalee2 scalee3" 2>&1 | tail -2

  # And the background parameters themselves, which is new: is the curvature
  # correlated with the scale parameters it is supposed to be protecting?
  if [[ $MODEL == expquad ]]; then
    cov "curv_vs_aem" "curv and A/e" \
      curv_BuJpsiKBkg_fit_bachelor_eta0 curv_BuJpsiKBkg_fit_bachelor_eta1 \
      scaleA0 scaleA1 scalee0 scalee1
  fi
done
echo "wrote to $OUT"
