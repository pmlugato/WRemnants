#!/usr/bin/env bash
# Combined fit: Z + J/psi + B+ -> J/psi K, one shared A/e/M basis.
#
# The whole point is the naming. merge_inputs.py correlates systematics BY NAME
# AND BY NOTHING ELSE, so the B+ tensor has to be written with
# --scaleParamNaming joint (A0..A23, e0..e23, M0..M23) rather than its
# standalone scaleA0.. names. With the standalone names the merge correlates
# ZERO parameters, the fit converges, and nothing downstream complains.
#
# The Z and J/psi tensors are those written by joint_mz_jpsi_fit.sh (the
# `repro` chain by default: Asimov, full 5D) and are read from <baseDir>/tensors.
#
# Usage:  joint_with_btojpsik.sh --baseDir DIR [options]
#   --baseDir DIR   base directory holding histograms/, tensors/, fitresults/
#                   and logs/ (required)
#   --stages LIST   comma-separated subset of 1 B+ tensor, 2 merge, 3 check,
#                   4 fit, 5 report, 6 plots                   [1,2,3,4,5,6]
#   --tag TAG       output label; every artefact carries it   [btojpsik_asimov]
#   --bhist FILE    B+ histmaker output, run at --scaleVariationEtaBins <nEta>
#                                 [<hists>/btojpsik_cvhnano_kpiA_kpi.hdf5]
#   --nEta N        shared nuisance eta basis                            [24]
#   --noJpsi        merge Z + B+ only, without the J/psi channels; the tag
#                   gets "_noJpsi"
#   --noB           merge Z + J/psi only (the nominal the others are compared
#                   to), with the same rabbit options; the tag gets "_noB"
#   --massBins N    B+ tensor mass bins over [5.05, 5.45]                [20]
#   --dataScale N   B+ Asimov yield x N, pinned yields follow            [1]
#   --mcStatScale F B+ template variances / F (more simulation)         [1]
#   --noZ           drop the Z channel (J/psi + B+, or J/psi alone with --noB)
#   --massShiftB N  B+ pole-mass POI, +/- N stored 5 MeV bins             [0]
#   --data          fit the observed data (blinded: only massShiftB* keeps rabbit's offset)
#   --jpsiModel M   J/psi param model: pou (per-cell, Asimov) | mu (signal strength, data)
#   --jpsiTensor T  J/psi tensor: repro (2 x 5D, Asimov) | full (1 channel, data-capable)
#   --zHist FILE    Z histmaker output behind the Z tensor, for the unit guard
#                                 [<hists>/mz_dilepton_..._joint_repro.hdf5]
#   --jpsiHist FILE J/psi histmaker output behind the J/psi tensor(s)
#                   [<hists>/dimuon_resonances_calinput_..._joint_<jpsiTensor>.hdf5]
#   --scaleVarParamScale "A=5 e=5"
#                   stretch the B+ A/e/M variations about the nominal, per
#                   parameter, at the tensor stage. Use with a 1/1/1 B+ histogram
#                   to reach the 5/5/1 unit of the full J/psi chain while the
#                   reweight stays linear                                [none]
#   --refTag TAG    tag of the fit the report compares against
#                                            [<tag without _noJpsi>_noB]
#   --plotDir DIR   base directory for the plot stage
#                       [~/public_html/mz/joint_mz_jpsi_bplus_<yymmdd>/<tag>]
#   --image PATH    singularity image
#
# The B+ arm is the primary standalone configuration: kpi histogram,
# peaking-b template, 24 free per-eta resolution nuisances at a 5 MeV
# reference, signal scale pinned to the data-fitted yield, 20 x 20 MeV over
# [5.05, 5.45]. Asimov is not blinded, so no --unblind is needed there.
#
# The joint fit of record:
#   joint_with_btojpsik.sh --baseDir DIR --tag res5_260914 --nEta 24 --massBins 20
# with --noJpsi / --noB for the Z + B+ and Z + J/psi references.
#
# Unit guard: stage 3 runs check_joint_tensor.py on
# the merge and stage 4 re-runs its provenance half and refuses to fit on a
# failure. One unit of A{i}/e{i}/M{i} must be the same shift in every channel.
# The repro Z/J/psi chain was built at scale 1/1/1 and the full one at 5/5/1, so
# the B+ histogram has to match whichever chain it is merged with.
#
# All configuration is by flag; TAG / BHIST / JOINT_N_ETA set in the environment
# are refused rather than ignored, so a stale export cannot change a fit.
set -uo pipefail

usage() { sed -n '2,/^set -uo/p' "${BASH_SOURCE[0]}" | sed '$d' | sed 's/^# \{0,1\}//'; }

for name in TAG BHIST JOINT_N_ETA; do
  if [[ -n ${!name+x} ]]; then
    echo "ERROR: $name is set in the environment; $(basename "$0") no longer reads it." >&2
    echo "Use --tag / --bhist / --nEta instead." >&2
    exit 2
  fi
done

WREM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BASE=""
IMAGE=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest

STAGES="1,2,3,4,5,6"
# The shared basis. 24, not 48: the J/psi histogram eta axis is 24 bins of 0.2,
# and at 48 adjacent nuisance pairs come out nearly degenerate (rho -0.81 / -0.75
# / -0.44 for A / e / M). The B+ histmaker must have been run with
# --scaleVariationEtaBins equal to this.
N_ETA=24
TAG="btojpsik_asimov"
BHIST=""
NOJPSI=0
NOB=0
NOZ=0
MASS_SHIFT_B=0
# --data: fit the observed histograms instead of the Asimov (-t 0 instead of -t -1).
# rabbit blinds every NOI on a data fit by a hashed offset; we release all of them
# except the B+ mass with --unblind, so only massShiftB* stays hidden.
DATA=0
UNBLIND_REGEX='^(?!massShiftB).*$'
# --jpsiModel pou|mu: 'pou' = per-cell AxisNormModel/AxisExpModel (the Asimov
# configuration; 9,552 free parameters that start far from real data and do not
# converge on it); 'mu' = rabbit's default per-channel signal strength, which is
# what the converged joint mZ + J/psi DATA fit used (joint_mz_jpsi_fit.sh, full mode).
JPSI_MODEL=pou
# --jpsiTensor repro|full: 'repro' = the two per-channel 5D tensors (230k bins each,
# need the 'pou' models, Asimov only); 'full' = joint_mz_jpsi_jpsi_full.hdf5, one
# channel, 5,760 bins, same 53.9M events and the same 72 A/e/M names -- the tensor
# the converged mZ + J/psi DATA fit used, with the 'mu' model.
JPSI_TENSOR=repro
MASS_BINS=20
REF_TAG=""
PLOTBASE=""
ZHIST=""
JPSIHIST=""
PARAM_SCALE=""
# Statistics scaling of the B+ Asimov, for projections only.
# --dataScale N multiplies the B+ Asimov yield by N; the pinned signal scale and
# its reference yield are multiplied by the same N so the model starts on the
# data. --mcStatScale F divides the B+ template variances by F (F = 1 is the
# simulation as produced; F = signal scale puts it on the data's footing).
DATA_SCALE=1
MC_STAT_SCALE=1
SIGNAL_SCALE_1FB=18.2433      # data-fitted signal / simulated signal, 2016 data
REF_YIELD_1FB=1327932         # the data-fitted signal yield behind it

need() { [[ $# -ge 2 ]] || { echo "$1 needs a value" >&2; exit 2; }; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --baseDir)  need "$@"; BASE=$2; shift 2 ;;
    --stages)   need "$@"; STAGES=$2; shift 2 ;;
    --tag)      need "$@"; TAG=$2; shift 2 ;;
    --bhist)    need "$@"; BHIST=$2; shift 2 ;;
    --nEta)     need "$@"; N_ETA=$2; shift 2 ;;
    --noJpsi)   NOJPSI=1; shift ;;
    --noB)      NOB=1; shift ;;
    --noZ)      NOZ=1; shift ;;
    --massShiftB)  need "$@"; MASS_SHIFT_B=$2; shift 2 ;;
    --data)     DATA=1; shift ;;
    --jpsiModel) need "$@"; JPSI_MODEL=$2; shift 2 ;;
    --jpsiTensor) need "$@"; JPSI_TENSOR=$2; shift 2 ;;
    --refTag)   need "$@"; REF_TAG=$2; shift 2 ;;
    --zHist)    need "$@"; ZHIST=$2; shift 2 ;;
    --scaleVarParamScale) need "$@"; PARAM_SCALE=$2; shift 2 ;;
    --jpsiHist) need "$@"; JPSIHIST=$2; shift 2 ;;
    --plotDir)  need "$@"; PLOTBASE=$2; shift 2 ;;
    --massBins) need "$@"; MASS_BINS=$2; shift 2 ;;
    --image)    need "$@"; IMAGE=$2; shift 2 ;;
    --dataScale)   need "$@"; DATA_SCALE=$2; shift 2 ;;
    --mcStatScale) need "$@"; MC_STAT_SCALE=$2; shift 2 ;;
    -h|--help)  usage; exit 0 ;;
    *) echo "unknown option '$1'" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n $BASE ]] || { echo "--baseDir <base directory> is required" >&2; exit 2; }
HISTS="$BASE/histograms"; TENSORS="$BASE/tensors"; FITS="$BASE/fitresults"; LOGS="$BASE/logs"
BHIST=${BHIST:-$HISTS/btojpsik_cvhnano_kpiA_kpi.hdf5}
SIGNAL_SCALE=$(python3 -c "print(f'{$SIGNAL_SCALE_1FB * $DATA_SCALE:.10g}')")
REF_YIELD=$(python3 -c "print(int(round($REF_YIELD_1FB * $DATA_SCALE)))")
if [[ $NOJPSI == 1 && $NOB == 1 ]]; then echo "--noJpsi and --noB together leave Z alone; not a configuration of this study" >&2; exit 2; fi
if [[ $NOZ == 1 && $NOJPSI == 1 ]]; then echo "--noZ and --noJpsi together leave B+ alone with a floating scale; not a configuration of this study" >&2; exit 2; fi
BASE_TAG=$TAG
[[ $NOJPSI == 1 ]] && TAG="${TAG}_noJpsi"
[[ $NOB == 1 ]] && TAG="${TAG}_noB"
[[ $NOZ == 1 ]] && TAG="${TAG}_noZ"
# The B+ tensor and the merge are shared between the Asimov and the data fit and
# do not depend on the J/psi param model, but they DO depend on which J/psi tensor
# is merged; so the merge name carries _jpsifull but not _data / _jpsimu.
MERGE_TAG=$TAG; [[ $JPSI_TENSOR == full ]] && MERGE_TAG="${MERGE_TAG}_jpsifull"
[[ $DATA == 1 ]] && TAG="${TAG}_data"
[[ $JPSI_MODEL == mu ]] && TAG="${TAG}_jpsimu"
[[ $JPSI_TENSOR == full ]] && TAG="${TAG}_jpsifull"
FIT_MODE="-t -1"; [[ $DATA == 1 ]] && FIT_MODE="-t 0 --unblind '$UNBLIND_REGEX'"
REF_TAG=${REF_TAG:-${BASE_TAG}_noB}
PLOTBASE=${PLOTBASE:-$HOME/public_html/mz/joint_mz_jpsi_bplus_$(date +%y%m%d)}
PLOTS="$PLOTBASE/$TAG"
run_stage() { [[ ",$STAGES," == *",$1,"* ]]; }

Z_TENSOR="$TENSORS/ZMassDilepton_nonTrigMuons_pt0_trigMuons_pt0_nonTrigMuons_eta0_trigMuons_eta0_mll_joint_repro/ZMassDilepton.hdf5"
JPSI_TENSORS="$TENSORS/joint_mz_jpsi_prompt_repro.hdf5 $TENSORS/joint_mz_jpsi_displaced_repro.hdf5"
[[ $JPSI_TENSOR == full ]] && JPSI_TENSORS="$TENSORS/joint_mz_jpsi_jpsi_full.hdf5"
case $JPSI_TENSOR in repro|full) ;; *) echo "--jpsiTensor must be repro or full" >&2; exit 2 ;; esac
[[ $NOJPSI == 1 ]] && JPSI_TENSORS=""
ZHIST=${ZHIST:-$HISTS/mz_dilepton_scetlib_dyturbo_CT18Z_N3p0LL_N2LO_Corr_joint_repro.hdf5}
JPSIHIST=${JPSIHIST:-$HISTS/dimuon_resonances_calinput_jpsi_2016PostVFP_etaBins_24_joint_${JPSI_TENSOR}.hdf5}
Z_PROJ="-m Project Z mll"
if [[ $NOZ == 1 ]]; then Z_TENSOR=""; Z_PROJ=""; fi
B_POSTFIX="joint${N_ETA}_${MERGE_TAG}"
B_TENSOR="$TENSORS/btojpsik_cvhnano_tensor_${B_POSTFIX}.hdf5"
MERGED="$TENSORS/joint_mz_jpsi_bplus_${MERGE_TAG}.hdf5"
FIT="$FITS/joint_mz_jpsi_bplus_${TAG}_fit.hdf5"

# B+ cell binning: 3 pT quantile bins x 24 eta (etaRebin 2, so 48/2) x 2 charge.
# 24 histogram eta bins matched 1:1 to the 24 nuisance groups -- a coarser
# histogram would leave several nuisances acting on the same bins.
B_ETA_REBIN=$((48 / N_ETA))
B_CELLS="fit_bachelor_pt,fit_bachelor_eta,fit_bachelor_charge"
B_MODELS="--paramModel AxisNormModel btojpsik_cvhnano BuJpsiK fit_bachelor_eta \
 --paramModel AxisNormModel btojpsik_cvhnano BuJpsiKMCPeak fit_bachelor_eta \
 --paramModel AxisExpPolyModel btojpsik_cvhnano BuJpsiKBkg jointCvh_mass $B_CELLS 2 fit_bachelor_eta"
B_PROJ="-m Project btojpsik_cvhnano jointCvh_mass"
if [[ $NOB == 1 ]]; then B_MODELS=""; B_PROJ=""; B_TENSOR=""; fi

# Inputs to the unit guard: one histogram per channel family in the merge.
GUARD_INPUTS=""
[[ -n $Z_TENSOR ]] && GUARD_INPUTS="$GUARD_INPUTS --zHist '$ZHIST' --zTensor '$Z_TENSOR'"
[[ -n $JPSI_TENSORS ]] && GUARD_INPUTS="$GUARD_INPUTS --jpsiHist '$JPSIHIST'"
[[ $NOB == 0 ]] && GUARD_INPUTS="$GUARD_INPUTS --bHist '$BHIST' --bTensor '$B_TENSOR'"
NFAMILIES=$(( (1-NOZ) + (1-NOJPSI) + (1-NOB) ))

# The J/psi param models and projections only exist when the J/psi tensors are
# in the merge; naming an absent channel is an error in rabbit, not a no-op.
JPSI_MODELS="--paramModel AxisNormModel JPsi_prompt sig_JPsi_prompt eta1,eta2,pt1,pt2 pou \
      --paramModel AxisExpModel JPsi_prompt background_JPsi_prompt mass eta1,eta2,pt1,pt2 pou \
      --paramModel AxisNormModel JPsi_displaced sig_JPsi_displaced eta1,eta2,pt1,pt2 pou \
      --paramModel AxisExpModel JPsi_displaced background_JPsi_displaced mass eta1,eta2,pt1,pt2 pou"
JPSI_PROJ="-m Project JPsi_prompt mass -m Project JPsi_displaced mass"
if [[ $NOJPSI == 1 ]]; then JPSI_MODELS=""; JPSI_PROJ=""; fi
if [[ $JPSI_MODEL == mu ]]; then JPSI_MODELS=""; fi
if [[ $JPSI_TENSOR == full ]]; then
  [[ $JPSI_MODEL == pou ]] && { echo "the full J/psi tensor has no per-cell axes for the pou models; use --jpsiModel mu" >&2; exit 2; }
fi   # channel names JPsi_prompt / JPsi_displaced are the same in both tensors, so JPSI_PROJ stands
case $JPSI_MODEL in pou|mu) ;; *) echo "--jpsiModel must be pou or mu" >&2; exit 2 ;; esac

run() {
  singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc \
    "cd '$WREM' && source setup.sh >/dev/null 2>&1 && $1"
}

mkdir -p "$TENSORS" "$FITS" "$LOGS"
echo "=== combined Asimov with B+, basis = $N_ETA eta bins, tag = $TAG$([[ $NOJPSI == 1 ]] && echo ', Z + B+ only') ==="
echo "    channels: Z=$((1-NOZ)) Jpsi=$((1-NOJPSI)) B+=$((1-NOB)); B+ mass POI shift bins: $MASS_SHIFT_B; $([[ $DATA == 1 ]] && echo "DATA fit, blinded except $UNBLIND_REGEX" || echo "Asimov")"
echo "    B+ statistics: data x $DATA_SCALE (signal scale $SIGNAL_SCALE, reference yield $REF_YIELD), MC-stat scale $MC_STAT_SCALE"

# ------------------------------------------------------------- 1. B+ tensor
if run_stage 1 && [[ $NOB == 0 ]]; then
  echo "--- [1/6] B+ tensor, joint naming"
  [[ -f "$BHIST" ]] || { echo "missing $BHIST -- run the histmaker at --scaleVariationEtaBins $N_ETA" >&2; exit 1; }
  run "/opt/venv/bin/python scripts/rabbit/btojpsik_cvhnano_tensor.py \
      -i '$BHIST' -o '$TENSORS' -p '$B_POSTFIX' \
      --ptBins 3 --etaRebin $B_ETA_REBIN \
      --scaleVariations --scaleVarScale 1.0 ${PARAM_SCALE:+--scaleVarParamScale $PARAM_SCALE} \
      --scaleParamNaming joint --jointEtaBins $N_ETA \
      --mcBkgTemplate peak --resolutionVariation 0.005 --resolutionAxes fit_bachelor_eta \
      --signalScale $SIGNAL_SCALE --signalScaleReferenceYield $REF_YIELD \
      --dataScale $DATA_SCALE --mcStatScale $MC_STAT_SCALE --massShiftB $MASS_SHIFT_B \
      --massLow 5.05 --massHigh 5.45 --massBins $MASS_BINS" 2>&1 | tee "$LOGS/joint_bplus_tensor_${TAG}.log"
  [[ ${PIPESTATUS[0]} == 0 ]] || { echo "B+ tensor stage FAILED -- stopping" >&2; exit 1; }
fi
[[ $NOB == 1 ]] && run_stage 1 && echo "--- [1/6] no B+ tensor in this configuration"

# ------------------------------------------------------------------ 2. merge
if run_stage 2; then
  echo "--- [2/6] merge"
  run "/opt/venv/bin/python rabbit/bin/merge_inputs.py \
      -t $Z_TENSOR $JPSI_TENSORS $B_TENSOR --sparse \
      --rename scripts/rabbit/lowmass_z_rename.py \
      -o '$TENSORS' --outname joint_mz_jpsi_bplus_${MERGE_TAG}" 2>&1 | tee "$LOGS/joint_bplus_merge_${MERGE_TAG}.log"
  [[ ${PIPESTATUS[0]} == 0 ]] || { echo "merge stage FAILED -- stopping" >&2; exit 1; }
fi

# ------------------------------------------------------------------ 3. check
if run_stage 3 && [[ $NOB == 0 ]] && [[ -n "$Z_TENSOR$JPSI_TENSORS" ]]; then
  echo "--- [3/6] how many parameters does the merge actually correlate?"
  run "/opt/venv/bin/python scripts/rabbit/check_bplus_correlation.py \
      --bplus '$B_TENSOR' --others $Z_TENSOR $JPSI_TENSORS \
      --merged '$MERGED' --rename scripts/rabbit/lowmass_z_rename.py \
      --nEta $N_ETA" 2>&1 | tee "$LOGS/joint_bplus_check_${TAG}.log"
fi
if run_stage 3 && [[ $NFAMILIES -ge 2 ]]; then
  echo "--- [3/6] unit guard: is one unit of A/e/M the same shift in every channel?"
  run "/opt/venv/bin/python scripts/rabbit/check_joint_tensor.py '$MERGED' \
      $GUARD_INPUTS --nEta $N_ETA" 2>&1 | tee "$LOGS/joint_bplus_unitguard_${TAG}.log"
  [[ ${PIPESTATUS[0]} == 0 ]] || { echo "unit guard FAILED -- not fitting this tensor" >&2; exit 1; }
fi

# -------------------------------------------------------------------- 4. fit
if run_stage 4; then
  if [[ $NFAMILIES -ge 2 ]]; then
    # Cheap (recorded arguments only), so it runs on every fit, including a
    # --stages 4 rerun that skipped stage 3.
    run "/opt/venv/bin/python scripts/rabbit/check_joint_tensor.py --provenanceOnly \
        $GUARD_INPUTS --nEta $N_ETA" 2>&1 | tee "$LOGS/joint_bplus_unitgate_${TAG}.log"
    [[ ${PIPESTATUS[0]} == 0 ]] || { echo "unit gate FAILED -- refusing to fit" >&2; exit 1; }
  fi
  echo "--- [4/6] combined fit"
  run "/opt/venv/bin/python rabbit/bin/rabbit_fit.py '$MERGED' \
      -o '$FITS' --outname joint_mz_jpsi_bplus_${TAG}_fit \
      $FIT_MODE --noChi2 --doImpacts --saveHists --saveHistsPerProcess \
      --computeHistErrors --hessianChunkSize 256 --verbose 4 \
      $JPSI_MODELS \
      $B_MODELS \
      $JPSI_PROJ \
      $Z_PROJ \
      $B_PROJ" 2>&1 | tee "$LOGS/joint_bplus_fit_${TAG}.log"
fi

# ----------------------------------------------------------------- 5. report
if run_stage 5 && [[ $NOB == 0 ]]; then
  echo "--- [5/6] report against $REF_TAG"
  REF="$FITS/joint_mz_jpsi_bplus_${REF_TAG}_fit.hdf5"
  run "/opt/venv/bin/python scripts/rabbit/report_bplus_gain.py '$FIT' --reference '$REF' \
      --nEta $N_ETA" 2>&1 | tee "$LOGS/joint_bplus_report_${TAG}.log"
fi
# ------------------------------------------------------------------ 6. plots
# In its own script so it can be re-run on an existing fit without editing
# this file while a chain is running (bash reads a running script
# incrementally, and an edit corrupts the running instance).
if run_stage 6; then
  echo "--- [6/6] plots -> $PLOTS"
  FLAGS=""; [[ $NOJPSI == 1 ]] && FLAGS="$FLAGS --noJpsi"; [[ $NOB == 1 ]] && FLAGS="$FLAGS --noB"
  [[ $DATA == 1 ]] && FLAGS="$FLAGS --data"
  if [[ $NOZ == 1 ]]; then
    FLAGS="$FLAGS --noZ"
    if [[ $MASS_SHIFT_B -gt 0 ]]; then FLAGS="$FLAGS --poi massShiftB$((MASS_SHIFT_B * 5))MeV --poiScale $((MASS_SHIFT_B * 5))"; else FLAGS="$FLAGS --poi ''"; fi
  fi
  "$WREM/scripts/rabbit/joint_bplus_plots.sh" "$FIT" "$PLOTS" "$N_ETA" "$TAG" $FLAGS
fi

echo "=== done ==="
