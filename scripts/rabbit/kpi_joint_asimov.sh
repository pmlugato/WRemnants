#!/usr/bin/env bash
# Section 6: the joint Z + J/psi + B+ Asimov, once per response-model
# configuration, against the Z + J/psi reference.
#
# Only the B+ arm changes between configurations. The Z and J/psi tensors are
# reused from the `repro` chain unchanged, which is what makes the comparison
# clean: any movement in sigma(mZ) is attributable to the B+ response model and
# to nothing else.
#
# Config C -- kpi on the J/psi legs as well -- is NOT here, and the reason is
# physics rather than effort. kpi returns the wrong SIGN above ~25 GeV, so it
# cannot back the Z arm at all; config C therefore reduces to "kpi on the J/psi
# legs only", which needs dimuon_resonances_calinput.py re-run with the
# --responseModel plumbing extended to it. Tracked as its own task.
#
# Usage: kpi_joint_asimov.sh [configs...]      default: all four
set -uo pipefail

CONFIGS=("$@")
[[ ${#CONFIGS[@]} -eq 0 ]] && CONFIGS=(muon_rec muon_cor kpi perleg)

WREM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CEPH=/ceph/submit/data/user/p/pmlugato/mz
HISTS="$CEPH/histograms"; FITS="$CEPH/fitresults"
REF="$FITS/joint_mz_jpsi_repro_fit.hdf5"

for cfg in "${CONFIGS[@]}"; do
  echo "================ joint Asimov, B+ arm = $cfg"
  date '+%H:%M'
  "$WREM/scripts/rabbit/joint_with_btojpsik.sh" \
      --bhist "$HISTS/btojpsik_cvhnano_kpiA_${cfg}.hdf5" --tag "kpi_${cfg}" --stages 1,2,3,4 2>&1 \
    | sed 's/\x1b\[[0-9;]*m//g' \
    | grep -iE "edmval|maximum number of iterations|of 72|FAIL|error|traceback|systematics" | head -10
done

echo "================ summary"
for cfg in "${CONFIGS[@]}"; do
  f="$FITS/joint_mz_jpsi_bplus_kpi_${cfg}_fit.hdf5"
  [[ -f $f ]] || { echo "$cfg: no fit"; continue; }
  echo "---- $cfg"
  "$WREM/scripts/rabbit/run_in_container.sh" \
    /opt/venv/bin/python "$WREM/scripts/rabbit/report_bplus_gain.py" \
      "$f" --reference "$REF" 2>&1 | tail -12
done
echo "JOINT KPI DONE"
