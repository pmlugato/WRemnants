#!/usr/bin/env bash
# rabbit_plot_cov.py readout of the A/e/M nuisance correlation, Asimov and data.
#
# One 12x12 matrix is unreadable with numbers printed on it, so the matrix is
# sliced three ways instead, each slice answering a different question:
#
#   eta{g}_aem   A/e/M inside one eta group  -> is A-e degenerate?
#   {p}_eta      one parameter across eta    -> do the eta groups share a mode?
#   Ae_all       A and e together, all eta   -> is the degeneracy local or global?
#
# Usage: cvhnano_correlation_plots.sh
set -uo pipefail

WREM_BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
MZ=${MZ:-/ceph/submit/data/user/p/pmlugato/mz}
OUTBASE=${OUTBASE:-$HOME/public_html/mz/cvhnano-aem/correlations}
IMAGE=${IMAGE:-/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest}

run() {
  singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc \
    "cd '$WREM_BASE' && source setup.sh >/dev/null 2>&1 && $1"
}

cov() {  # cov <fit.hdf5> <outdir> <tag> <xlabel> <params...>
  local fit=$1 out=$2 tag=$3 lab=$4; shift 4
  mkdir -p "$out"
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$fit' -v 4 -o '$out' \
       --correlation --params $* --showNumbers -p '$tag' \
       --xlabel '$lab' --ylabel '$lab'" 2>&1 | tail -3
}

declare -A FITS=(
  [asimov]="$MZ/fitresults/cvhnano_full_exp_asimov_eta12_aem_full_fit.hdf5"
  [data]="$MZ/fitresults/cvhnano_full_exp_data_eta12_normgrouped_aem_full_fit.hdf5"
)

for which in asimov data; do
  FIT=${FITS[$which]}
  OUT="$OUTBASE/$which"
  [[ -f $FIT ]] || { echo "missing $FIT" >&2; continue; }
  echo "=== $which : $FIT"

  # A/e/M within one eta group: the degeneracy the 2D scan shows, per group.
  for g in 0 1 2 3; do
    cov "$FIT" "$OUT" "eta${g}_aem" "A/e/M  eta group $g" \
        "scaleA$g" "scalee$g" "scaleM$g"
  done

  # One parameter across eta groups.
  for p in A e M; do
    cov "$FIT" "$OUT" "${p}_across_eta" "$p by eta group" \
        "scale${p}0" "scale${p}1" "scale${p}2" "scale${p}3"
  done

  # A and e together, every group: is the degeneracy block-diagonal? No printed
  # numbers -- at 8x8 they overlap into an unreadable smear, and the question
  # here is a pattern (does A_i correlate only with e_i) which colour answers.
  mkdir -p "$OUT"
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -v 4 -o '$OUT' \
       --correlation -p 'Ae_all' --xlabel 'A (0-3) then e (0-3)' \
       --ylabel 'A (0-3) then e (0-3)' \
       --params scaleA0 scaleA1 scaleA2 scaleA3 scalee0 scalee1 scalee2 scalee3" 2>&1 | tail -3

  # The full block, no numbers: the question there is "is it diagonal", and
  # colour answers it where 144 printed values do not.
  mkdir -p "$OUT"
  run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -v 4 -o '$OUT' \
       --correlation -p 'aem_all' --xlabel 'A/e/M, all eta groups' \
       --ylabel 'A/e/M, all eta groups' \
       --params scaleA0 scaleA1 scaleA2 scaleA3 scalee0 scalee1 scalee2 scalee3 \
                scaleM0 scaleM1 scaleM2 scaleM3" 2>&1 | tail -3
done

echo; echo "plots under $OUTBASE"
