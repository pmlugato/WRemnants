#!/usr/bin/env bash
# Prefit, postfit, uncertainty and correlation readout for every CONVERGED
# Asimov configuration of the staged model build.
#
# Asimov only, and that is not a limitation of this script: no data fit in this
# analysis converges, so a data correlation matrix would be a Hessian at a
# non-stationary point, which is not a covariance.
set -uo pipefail
WREM_BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
MZ=${MZ:-/ceph/submit/data/user/p/pmlugato/mz}
OUT=${OUT:-$HOME/public_html/mz/cvhnano-aem/stages}
IMAGE=${IMAGE:-/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest}
CH=btojpsik_cvhnano; MASS=jointCvh_mass

run() { singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc \
    "cd '$WREM_BASE' && source setup.sh >/dev/null 2>&1 && $1" 2>&1 | tail -2; }

CONFIGS=${CONFIGS:-"n2_full n2_mccomb_full n2_mcall_full n2_mcsplit_full n2_res00089_full n2_res00089c_full"}
for TAG in $CONFIGS; do
  BASE="cvhnano_full_exppoly_asimov_eta12_normgrouped_aem_${TAG}"
  FIT="$MZ/fitresults/${BASE}_fit.hdf5"; PR="$MZ/fitresults/${BASE}_plotready.hdf5"
  [[ -f $FIT && -f $PR ]] || { echo "skip $TAG (missing)"; continue; }
  D="$OUT/$TAG"; mkdir -p "$D"
  echo "=== $TAG"
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists.py '$PR' -o '$D' -p prefit --prefit -m Project $CH $MASS --dataName Asimov"
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists.py '$PR' -o '$D' -p postfit -m Project $CH $MASS --dataName Asimov"
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists_uncertainties.py '$PR' -o '$D' -p unc -m Project $CH $MASS"
  for g in 0 1 2 3; do
    run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -o '$D' --correlation \
         --params scaleA${g} scalee${g} scaleM${g} --showNumbers -p eta${g}_aem"
  done
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -o '$D' --correlation -p Ae_all \
       --params scaleA0 scaleA1 scaleA2 scaleA3 scalee0 scalee1 scalee2 scalee3"
done
echo "wrote to $OUT"
