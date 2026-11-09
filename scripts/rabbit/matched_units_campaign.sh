#!/usr/bin/env bash
# Refits with the B+ A/e/M unit matched to the Z and J/psi (notes/aem_unit_mismatch.md, step 2).
#
# The two Z/J/psi chains carry different units (repro 1/1/1, full 5/5/1). The B+
# is built once, at 1/1/1, where the kaon reweight is linear:
#   kpiC_fitscale1  --fitMuonScaleAndResolution, scale 1/1/1
#     repro Z + repro J/psi : used as is
#     full J/psi (mB fits)  : A and e stretched x5 at the tensor stage
#                             (--scaleVarParamScale "A=5 e=5"), M x1
# kpiC_fitscale5 (built directly at 5/5/1) is NOT used for fits: its e response
# saturates at low kaon pT (aem_unit_mismatch.md 9.3); it is the linearity test.
# Every merged fit goes through joint_with_btojpsik.sh, whose unit guard refuses
# a mismatched merge. Each step mirrors the configuration of the fit it is
# compared with, changing only the B+ histogram.
#
# Usage: matched_units_campaign.sh STEP [STEP ...]
#   bplus      standalone B+ Asimov, 1/1/1 tensor         (vs alllegs_bplus_asimov)
#   mz         Z + J/psi + B+ Asimov, repro chain        (vs res5_260914, ref _noB)
#   mBasimov   J/psi + B+, mB floating, Asimov           (vs realmc_C_asimov)
#   mBdata     the same on DATA, mB blinded              (vs realmc_A_data)
#              Run only once the Asimov fits make sense.
set -uo pipefail
WREM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$WREM"
CEPH=/ceph/submit/data/user/p/pmlugato/mz
H=$CEPH/histograms; T=$CEPH/tensors; F=$CEPH/fitresults; L=$CEPH/logs
IMAGE=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest
B1=$H/btojpsik_cvhnano_kpiC_fitscale1.hdf5
DRV=scripts/rabbit/joint_with_btojpsik.sh
B_MODELS="--paramModel AxisNormModel btojpsik_cvhnano BuJpsiK fit_bachelor_eta \
 --paramModel AxisNormModel btojpsik_cvhnano BuJpsiKMCPeak fit_bachelor_eta \
 --paramModel AxisExpPolyModel btojpsik_cvhnano BuJpsiKBkg jointCvh_mass fit_bachelor_pt,fit_bachelor_eta,fit_bachelor_charge 2 fit_bachelor_eta"

run() {
  singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc \
    "cd '$WREM' && source setup.sh >/dev/null 2>&1 && $1"
}
[[ $# -ge 1 ]] || { sed -n '2,/^set -uo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit 2; }
for f in "$B1"; do [[ -r $f ]] || { echo "missing $f" >&2; exit 1; }; done

for step in "$@"; do
  case $step in
    mz)
      # Stage 1 writes the B+ tensor btojpsik_cvhnano_tensor_joint24_matched_mz.hdf5,
      # which the standalone step reuses.
      $DRV --tag matched_mz --bhist "$B1" --refTag res5_260914_noB --stages 1,2,3,4,5 \
        2>&1 | tee "$L/matched_mz.log" ;;
    mBasimov)
      $DRV --tag matched_mB --noZ --massShiftB 1 --jpsiModel mu --jpsiTensor full \
        --bhist "$B1" --scaleVarParamScale "A=5 e=5" --stages 1,2,3,4 2>&1 | tee "$L/matched_mB_asimov.log" ;;
    mBdata)
      # Blinded: the driver's --data runs -t 0 --unblind '^(?!massShiftB).*$', so
      # the B+ mass keeps rabbit's hashed offset. Reuses the Asimov chain's merge.
      $DRV --tag matched_mB --noZ --massShiftB 1 --jpsiModel mu --jpsiTensor full \
        --bhist "$B1" --scaleVarParamScale "A=5 e=5" --data --stages 4 2>&1 | tee "$L/matched_mB_data.log" ;;
    bplus)
      # Only the 1/1/1 tensor: the mB-chain tensor carries the floating B mass,
      # and B+ alone with both mB and the scale free is degenerate.
      for pair in "1:$T/btojpsik_cvhnano_tensor_joint24_matched_mz.hdf5"; do
        sc=${pair%%:*}; ten=${pair#*:}
        [[ -r $ten ]] || { echo "missing $ten -- run the mz / mBasimov step first" >&2; exit 1; }
        run "/opt/venv/bin/python rabbit/bin/rabbit_fit.py '$ten' -o '$F' \
            --outname matched${sc}_bplus_asimov_fit -t -1 --noChi2 --doImpacts --saveHists \
            --computeHistErrors --verbose 4 $B_MODELS -m Project btojpsik_cvhnano jointCvh_mass" \
          2>&1 | tee "$L/matched${sc}_bplus_asimov.log"
      done ;;
    *) echo "unknown step '$step'" >&2; exit 2 ;;
  esac
done
