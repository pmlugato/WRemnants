#!/bin/bash
# Joint mZ + J/psi fit: histmakers -> tensors -> merge -> check -> fit -> plots.
#
# Correlation between the Z and J/psi channels happens by parameter name and by
# nothing else (see scripts/rabbit/lowmass_z_rename.py). The settings below are
# not defaults and must not be dropped: every one of them is a shared assumption
# whose mismatch produces a wrong mZ without producing an error.
# check_joint_tensor.py verifies them from the recorded arguments before the fit
# is allowed to run.
#
#   --etaBins <nEta>               same on BOTH histmakers (--nEta, default
#                                  24), and the --nEta the check is run with
#   --scale_A 5.0 --scale_e 5.0    same on BOTH histmakers. The J/psi histmaker
#                                  overrides BOTH to 5.0 via set_parser_default;
#                                  mz_dilepton inherits 1.0 for both. Missing
#                                  --scale_A was caught by check_joint_tensor.py
#                                  after being missed on a source read.
#   --fitMuonScaleAndResolution    on BOTH histmakers
#   --calibrationStatScaling 1.0   setupRabbit default is 2.1, Z side only
#   --noClosureSysts               the closure systs cover a correction that is
#                                  no longer applied
#
# Usage:
#   joint_mz_jpsi_fit.sh repro            # havyn's exact recipe, full 5D, sparse
#   joint_mz_jpsi_fit.sh pilot            # 1-file Z, full J/psi, coarse binning
#   joint_mz_jpsi_fit.sh full             # everything
#   joint_mz_jpsi_fit.sh pilot 4,5,6      # only the listed stages
#
# Options, after the positionals:
#   --baseDir DIR   base directory for histograms/, tensors/, fitresults/ and
#                   logs/ (required)
#   --nEta N        nuisance eta basis for both histmakers and the check  [24]
#   --nThreads N    histmaker threads (-j)                               [160]
#
# Stages: 1 Z histmaker, 2 J/psi histmaker, 3 Z tensor, 4 J/psi tensor,
#         5 merge + check, 6 fit, 7 plots.

set -euo pipefail

MODE=pilot
STAGES=1,2,3,4,5,6,7
if [[ $# -gt 0 && $1 != --* ]]; then MODE=$1; shift; fi
if [[ $# -gt 0 && $1 != --* ]]; then STAGES=$1; shift; fi

case "$MODE" in
  pilot|full|repro) ;;
  *) echo "usage: $0 {pilot|full|repro} [stages] --baseDir DIR [--nEta N] [--nThreads N]" >&2; exit 2 ;;
esac

# Configuration is by flag only; a stale export must not change a fit silently.
for name in JOINT_N_ETA NTHREADS; do
  if [[ -n ${!name+x} ]]; then
    echo "ERROR: $name is set in the environment; $(basename "$0") does not read it." >&2
    echo "Use --nEta / --nThreads instead." >&2
    exit 2
  fi
done

# Nuisance eta binning, shared by BOTH histmakers. The calibration file is
# native 48 bins of 0.1, but the J/psi histogram eta axis is 24 bins of 0.2, so
# at 48 each histogram bin spans two nuisance bins and adjacent pairs come out
# nearly degenerate (measured rho = -0.81 / -0.75 / -0.44 for A / e / M). 24
# aligns the basis 1:1 with what the data resolves.
BASE=""
N_ETA=24
# The histmakers are I/O-bound on the NanoAOD reads, so more threads than
# cores-in-use still helps; lower it on a busy node.
NTHREADS=160
while [[ $# -gt 0 ]]; do
  case "$1" in
    --baseDir)  [[ $# -ge 2 ]] || { echo "--baseDir needs a value" >&2; exit 2; }; BASE=$2; shift 2 ;;
    --nEta)     [[ $# -ge 2 ]] || { echo "--nEta needs a value" >&2; exit 2; }; N_ETA=$2; shift 2 ;;
    --nThreads) [[ $# -ge 2 ]] || { echo "--nThreads needs a value" >&2; exit 2; }; NTHREADS=$2; shift 2 ;;
    *) echo "unknown option '$1'" >&2; exit 2 ;;
  esac
done

run_stage() { [[ ",$STAGES," == *",$1,"* ]]; }

WREM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
[[ -n $BASE ]] || { echo "--baseDir <base directory> is required" >&2; exit 2; }
HISTS="$BASE/histograms"
TENSORS="$BASE/tensors"
FITS="$BASE/fitresults"
LOGS="$BASE/logs"
PLOTS="$HOME/public_html/mz"
TAG="joint_${MODE}"

mkdir -p "$HISTS" "$TENSORS" "$FITS" "$LOGS"
mkdir -p "$PLOTS/joint_mz_jpsi_prefit_postfit" \
         "$PLOTS/joint_mz_jpsi_pulls_impacts" \
         "$PLOTS/joint_mz_jpsi_covariance" \
         "$PLOTS/joint_mz_jpsi_datamc"

# Everything runs inside the container with setup.sh sourced.
#
# Pass the command as ARGUMENTS via run_in_container.sh, never as one string to
# `bash -c`. The container's runscript rebuilds the command line by wrapping
# each argument in double quotes and eval-ing the result, so a command handed
# over as a single string is word-split a second time inside the container.
# Verified: with this form, --unblind 'sig_JPsi.*' arrives as one argument and
# --postfix 'x y' keeps its space.
#
# The one rule that remains: no argument may contain a literal double quote,
# since that closes the runscript's own quoting. Use single quotes inside any
# python -c snippet.
sing() {
  "$WREM/scripts/ci/run_with_singularity.sh" "$WREM/scripts/rabbit/run_in_container.sh" "$@"
}

# ---------------------------------------------------------------------------
# repro: havyn's recipe, verbatim. Keeps the full 5D
# binning on BOTH channels and controls size with --sparse plus an MC-event
# threshold, rather than projecting p_T away. That matters physically: the
# correction is dk/k = A - e*k with k = 1/pT, so A and e are told apart ONLY by
# their p_T dependence. Projecting p_T out drove rho(A,e) to 0.98 against
# havyn's 0.90 and inflated the calibration impact.
if [[ "$MODE" == repro ]]; then
  Z_MAXFILES=""
  Z_AXES="nonTrigMuons_pt0 trigMuons_pt0 nonTrigMuons_eta0 trigMuons_eta0 mll"
  Z_BINNING="--eta 24 -2.4 2.4 --pt 3 26. 60."
  Z_FITVAR="nonTrigMuons_pt0-trigMuons_pt0-nonTrigMuons_eta0-trigMuons_eta0-mll"
  Z_EXTRA="--noi zmass --massConstraintModeZ unconstrained"
  # No --realData: havyn's slide is "Asimov first results".
  Z_DATA=""
  JPSI_EXTRA="--pixelMultiplicityStat"
  JPSI_BINNING="--sparse --mcEventThresh 500"
  MERGE_EXTRA="--sparse"
  FIT_TOYS="-1"          # Asimov
else
  Z_AXES="etaAbsEta mll"
  Z_BINNING=""
  Z_FITVAR="etaAbsEta-mll"
  Z_EXTRA=""
  Z_DATA="--realData"
  JPSI_EXTRA=""
  # NOT --massOnly: collapsing eta makes the scale nuisances degenerate. NOT the
  # full 5D either: 46080 bins/channel puts nbins*nsyst*nproc past int32.
  JPSI_BINNING="--massRebin 5 --projectAxes eta1 eta2"
  MERGE_EXTRA=""
  FIT_TOYS="0"
  if [[ "$MODE" == pilot ]]; then Z_MAXFILES="--maxFiles 1"; else Z_MAXFILES=""; fi
fi

# Both histmakers at 1.0, matching havyn. The J/psi histmaker's own defaults are
# 5.0 for BOTH scale_A and scale_e (set_parser_default) while
# mz_dilepton inherits 1.0, so scale_A/scale_e must be passed explicitly to the
# J/psi side to bring it down.
SCALES="--scale_A 1.0 --scale_e 1.0"

echo "=== joint mZ + J/psi chain, mode=$MODE stages=$STAGES ==="

# ---------------------------------------------------------------- 1. Z hists
if run_stage 1; then
  echo "--- [1/7] mz_dilepton"
  sing python $WREM/scripts/histmakers/mz_dilepton.py \
      --axes $Z_AXES $Z_BINNING $Z_MAXFILES -j $NTHREADS \
      --fitMuonScaleAndResolution $SCALES --etaBins $N_ETA \
      -o $HISTS --postfix $TAG 2>&1 | tee "$LOGS/joint_mz_jpsi_z_${MODE}.log"
fi
Z_HIST=$(ls -t "$HISTS"/mz_dilepton_*_"$TAG".hdf5 2>/dev/null | head -1 || true)

# ------------------------------------------------------------ 2. J/psi hists
if run_stage 2; then
  echo "--- [2/7] dimuon_resonances_calinput"
  sing python $WREM/scripts/histmakers/dimuon_resonances_calinput.py \
      --resonance jpsi --triggers dimuon20_jpsi doublemu4_jpsitrk_displaced \
      --fitMuonScaleAndResolution $SCALES --etaBins $N_ETA $JPSI_EXTRA \
      -j $NTHREADS -o $HISTS --postfix $TAG 2>&1 | tee "$LOGS/joint_mz_jpsi_calinput_${MODE}.log"
fi
JPSI_HIST=$(ls -t "$HISTS"/dimuon_resonances_calinput_*_"$TAG".hdf5 2>/dev/null | head -1 || true)

# -------------------------------------------------------------- 3. Z tensor
Z_TENSOR_DIR="$TENSORS/ZMassDilepton_$(echo $Z_FITVAR | tr - _)_$TAG"
if run_stage 3; then
  echo "--- [3/7] setupRabbit"
  [[ -n "$Z_HIST" ]] || { echo "no Z histmaker output found" >&2; exit 1; }
  sing python $WREM/scripts/rabbit/setupRabbit.py -i $Z_HIST \
      --fitvar $Z_FITVAR $Z_DATA $Z_EXTRA \
      --noClosureSysts --calibrationStatScaling 1.0 \
      -o $TENSORS --postfix $TAG 2>&1 | tee "$LOGS/joint_mz_jpsi_ztensor_${MODE}.log"
fi
Z_TENSOR="$Z_TENSOR_DIR/ZMassDilepton.hdf5"

# ---------------------------------------------------------- 4. J/psi tensor(s)
# repro writes ONE TENSOR PER CHANNEL, as havyn does. The channels stay separate
# through the merge either way; splitting them here keeps each sparse tensor's
# event threshold applied to its own cells.
if [[ "$MODE" == repro ]]; then
  JPSI_TENSORS="$TENSORS/joint_mz_jpsi_prompt_$MODE.hdf5 $TENSORS/joint_mz_jpsi_displaced_$MODE.hdf5"
else
  JPSI_TENSORS="$TENSORS/joint_mz_jpsi_jpsi_$MODE.hdf5"
fi
if run_stage 4; then
  echo "--- [4/7] lowmass_tensor"
  [[ -n "$JPSI_HIST" ]] || { echo "no J/psi histmaker output found" >&2; exit 1; }
  if [[ "$MODE" == repro ]]; then
    for CH in prompt displaced; do
      sing python $WREM/scripts/rabbit/lowmass_tensor.py \
          --calinputHdf5 $JPSI_HIST $JPSI_BINNING \
          --jpsiChannels JPsi_$CH \
          -o $TENSORS --outname joint_mz_jpsi_${CH}_$MODE \
          2>&1 | tee -a "$LOGS/joint_mz_jpsi_jtensor_${MODE}.log"
    done
  else
    sing python $WREM/scripts/rabbit/lowmass_tensor.py \
        --calinputHdf5 $JPSI_HIST $JPSI_BINNING \
        -o $TENSORS --outname joint_mz_jpsi_jpsi_$MODE 2>&1 | tee "$LOGS/joint_mz_jpsi_jtensor_${MODE}.log"
  fi
fi

# ------------------------------------------------------- 5. merge and check
MERGED="$TENSORS/joint_mz_jpsi_merged_$MODE.hdf5"
if run_stage 5; then
  echo "--- [5/7] merge"
  sing python $WREM/rabbit/bin/merge_inputs.py \
      -t $Z_TENSOR $JPSI_TENSORS $MERGE_EXTRA \
      --rename $WREM/scripts/rabbit/lowmass_z_rename.py \
      -o $TENSORS --outname joint_mz_jpsi_merged_$MODE 2>&1 | tee "$LOGS/joint_mz_jpsi_merge_${MODE}.log"

  echo "--- [5/7] check the correlation contract"
  # Gates the fit. A failure here means the tensor correlates something other
  # than what it claims to, and the resulting mZ would be wrong-but-plausible.
  sing python $WREM/scripts/rabbit/check_joint_tensor.py $MERGED \
      --zHist $Z_HIST --jpsiHist $JPSI_HIST --zTensor $Z_TENSOR \
      --nEta $N_ETA --zChannel Z \
      2>&1 | tee "$LOGS/joint_mz_jpsi_check_${MODE}.log"

  echo "--- [5/7] debug_inputdata"
  sing python $WREM/rabbit/bin/rabbit_debug_inputdata.py $MERGED \
      2>&1 | tee "$LOGS/joint_mz_jpsi_debug_${MODE}.log" || true
fi

# ------------------------------------------------------------------ 6. fit
FITRESULT="$FITS/joint_mz_jpsi_${MODE}_fit.hdf5"
if run_stage 6; then
  echo "--- [6/7] fit"
  # --noChi2 is required, not optional. save_hists computes a chi2 whose
  # covariance is pred_cov + data_cov, and lowmass_tensor.py zeroes structurally
  # empty cells by design, so those rows are identically zero in both terms and
  # tf.linalg.solve raises "Input matrix is not invertible" -- in the PREFIT
  # hist save, before the minimiser ever runs. havyn passes it too.
  if [[ "$MODE" == repro ]]; then
    # havyn's command verbatim. The J/psi signal gets a per-cell floating
    # normalisation and the background an exponential in mass, both as `pou`
    # (parameters of uninterest) rather than POIs -- so the only POI/NOI in the
    # fit is the Z mass. The Project mappings are what the plots are made from;
    # there is no BaseMapping here.
    sing python $WREM/rabbit/bin/rabbit_fit.py $MERGED \
        -o $FITS --outname joint_mz_jpsi_${MODE}_fit \
        -t -1 --noChi2 --doImpacts --saveHists --saveHistsPerProcess \
        --computeHistErrors --hessianChunkSize 256 --verbose 4 \
        --paramModel AxisNormModel JPsi_prompt sig_JPsi_prompt eta1,eta2,pt1,pt2 pou \
        --paramModel AxisExpModel JPsi_prompt background_JPsi_prompt mass eta1,eta2,pt1,pt2 pou \
        --paramModel AxisNormModel JPsi_displaced sig_JPsi_displaced eta1,eta2,pt1,pt2 pou \
        --paramModel AxisExpModel JPsi_displaced background_JPsi_displaced mass eta1,eta2,pt1,pt2 pou \
        -m Project JPsi_prompt mass -m Project JPsi_prompt eta1 mass \
        -m Project JPsi_displaced mass -m Project JPsi_displaced eta1 mass \
        -m Project Z mll -m Project Z nonTrigMuons_eta0 mll \
        2>&1 | tee "$LOGS/joint_mz_jpsi_fit_${MODE}.log"
  else
    # --unblind on the J/psi signal strengths only: they are POIs under the
    # default Mu model but are not the measurement, and their multiplicative
    # exp(N(0,5)) blinding offset would start the fit orders of magnitude from
    # the data. mZ stays blinded, which costs nothing: the offset is additive on
    # the NOI, so uncertainty, covariance, pulls and constraints are unaffected.
    sing python $WREM/rabbit/bin/rabbit_fit.py $MERGED \
        -t $FIT_TOYS --doImpacts --saveHists --saveHistsPerProcess --computeHistErrors \
        --computeVariations --noChi2 --unblind 'sig_JPsi.*' \
        -m BaseMapping \
        -o $FITS --outname joint_mz_jpsi_${MODE}_fit 2>&1 | tee "$LOGS/joint_mz_jpsi_fit_${MODE}.log"
  fi
fi

# ---------------------------------------------------------------- 7. plots
if run_stage 7; then
  echo "--- [7/7] plots"
  # -m BaseMapping is REQUIRED: rabbit_plot_hists loops over `args.mapping`,
  # which defaults to empty, so without it the script exits 0 having plotted
  # nothing at all. A silent success is the worst failure mode here.
  if [[ "$MODE" == repro ]]; then
    MAPS="-m Project"
  else
    MAPS="-m BaseMapping"
  fi
  for FT in "" "--prefit"; do
    sing python $WREM/rabbit/bin/rabbit_plot_hists.py $FITRESULT \
        -o $PLOTS/joint_mz_jpsi_prefit_postfit $MAPS $FT \
        2>&1 | tee -a "$LOGS/joint_mz_jpsi_plots_${MODE}.log" || true
  done
  sing python $WREM/rabbit/bin/rabbit_plot_pulls_and_impacts.py $FITRESULT \
      -o $PLOTS/joint_mz_jpsi_pulls_impacts 2>&1 | tee -a "$LOGS/joint_mz_jpsi_plots_${MODE}.log" || true
  sing python $WREM/rabbit/bin/rabbit_print_impacts.py $FITRESULT \
      2>&1 | tee "$LOGS/joint_mz_jpsi_impacts_${MODE}.txt" || true
  # --params takes exact names, not regexes; a central-eta block stays readable.
  CLO=$(( N_ETA / 2 - 4 )); CHI=$(( N_ETA / 2 + 3 ))
  COVPARAMS="massShiftZ100MeV $(for i in $(seq $CLO $CHI); do echo -n "A$i e$i M$i "; done)"
  sing python $WREM/rabbit/bin/rabbit_plot_cov.py $FITRESULT \
      -o $PLOTS/joint_mz_jpsi_covariance --params $COVPARAMS \
      2>&1 | tee -a "$LOGS/joint_mz_jpsi_plots_${MODE}.log" || true
  # Joint-fit-specific views the central rabbit scripts do not produce:
  # constraint vs eta, mZ correlation vs eta, adjacent-pair correlation, and
  # the grouped impacts as a chart.
  sing python $WREM/scripts/plotting/joint_mz_jpsi_summary.py $FITRESULT \
      -o $PLOTS/joint_mz_jpsi_pulls_impacts --tag $MODE \
      2>&1 | tee -a "$LOGS/joint_mz_jpsi_plots_${MODE}.log" || true
fi

echo "=== done ==="
echo "  Z hist     : $Z_HIST"
echo "  J/psi hist : $JPSI_HIST"
echo "  merged     : $MERGED"
echo "  fit        : $FITRESULT"
echo "  plots      : $PLOTS/joint_mz_jpsi_*"
