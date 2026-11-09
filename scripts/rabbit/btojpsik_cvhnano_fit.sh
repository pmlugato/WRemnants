#!/usr/bin/env bash
# Tensor -> fit -> plot-ready chain for the CVH NanoAOD B+ -> J/psi K channel.
#
# The full cell binning is the primary configuration. The coarse ones are
# diagnostics for localising a failure, not a ladder to climb first; the one
# exception is `inclusive`, which is a deliverable in its own right because the
# signal yield in data comes from it.
#
# A fit that does not converge is an acceptable outcome. This script therefore
# does not stop on a failed fit: it records the failure and exits 0, so the
# attribution step has every configuration's output to work from.
#
# Usage:
#   btojpsik_cvhnano_fit.sh <histmaker.hdf5> <config> [options] [-- extra tensor args...]
#
#   config    full | ptcharge | inclusive
#
# Everything is a command-line flag. This script used to read ~25 knobs from the
# environment (BKGMODEL=..., RESVAR=..., ...); that was retired on 2026-09-11
# because a stale export in the shell silently changed a fit with nothing on the
# command line to show it. If any of the old names is set, the script refuses to
# run and says which flag replaces it.
#
# Options (defaults in brackets):
#   --bkgModel M         exp | expquad | exppoly | bernstein            [exp]
#   --polyOrder N        order for --bkgModel exppoly. Order 1 is exp, order 2
#                        is expquad; measured gain is 2.00x at order 2, 2.03x at
#                        3, 2.07x at 4, nothing at 5.                     [2]
#   --observed           fit the observation instead of an Asimov dataset. Runs
#                        UNBLINDED: rabbit multiplies every normalisation by
#                        exp(N(0, 5)) on an observed fit unless --unblind is
#                        passed, which put every B+ data fit before 2026-09-14
#                        up to e^10 off at the start (prefit loss 986k blinded,
#                        46k unblinded) and made trust-krylov stall. A
#                        calibration channel has nothing to blind. Applied to
#                        the fit AND the plot-ready pass.
#   --scaleVar           add the A/e/M nuisances from the histmaker's variation
#                        tensor (off reproduces the old no-systematics fit)
#   --scaleVarScale X    inflate every variation about the nominal      [1.0]
#   --skipScaleParams S  space-separated A/e/M names to leave out       [none]
#   --scan2D "pA pB"     run a 2D likelihood scan on this pair          [none]
#   --etaRebin N         integer rebin of the 48-bin eta axis             [1]
#   --ptBins N           number of pT quantile bins                       [3]
#   --slopeAxes CSV      cell axes the exponential slope varies over
#                                                            [fit_bachelor_eta]
#   --curvAxes CSV       cell axes the expquad curvature varies over. Default
#                        follows --slopeAxes; pass an empty string ("") for one
#                        global curvature.
#   --normAxes CSV       cell axes the signal normalisation varies over
#                        (empty = every cell axis)                       [all]
#   --mcBkg M            none | all | comb | split | peak -- add the simulation's
#                        own background shape as template process(es) ALONGSIDE
#                        the analytic one. The analytic component is always kept:
#                        the simulation has no prompt and no non-b fake
#                        background, so it cannot describe the data by itself.
#                        `peak` is the peaking-b template (gen categories 1-5)
#                        alone, without the noisy combinatorial template. [none]
#   --mcBkgNormAxes CSV  cell axes the MC template normalisations vary over
#                        (default: same as --normAxes)
#   --bkgSeed M          flat | expSideband -- shape of the analytic background
#                        template, i.e. the optimiser's start; the fitted model
#                        class is identical (a linear exponent is absorbed by
#                        the poly's first order).                        [flat]
#   --resVar W           signal resolution smearing width in GeV; 0 disables.
#                        Use ~0.003 (matched-purity fit, openspec 17.5); NOT
#                        the 0.0089 this once said.                       [0]
#   --resConstrained     give the resolution nuisance a unit Gaussian prior
#                        (default: unconstrained, the data determines it)
#   --resAxes CSV        cell axes the resolution nuisance varies over. Empty
#                        keeps one global sigRes; fit_bachelor_eta gives one free
#                        nuisance per eta bin (sigRes_eta{i}, the physics choice);
#                        fit_bachelor_eta,fit_bachelor_pt one per (eta, pT).
#                        Needs --resVar.                                   [""]
#   --signalScale X      pin the global signal template scale instead of the
#                        capped seed estimate. Set it to (data-fitted signal
#                        yield)/(unscaled MC signal in [5.05, 5.45]) so the
#                        Asimov carries the data's signal yield; the seed cap
#                        put every earlier Asimov 5-20% low.           [seed]
#   --refYield Y         metadata: the yield --signalScale was derived from
#   --massLow X --massHigh Y --massBins N
#                        fit window and bin count (writer defaults: 5.20, 5.35,
#                        10). The primary of add-cvhnano-per-eta-resolution is
#                        5.05 5.45 20; its fallback is 5.20 5.40 10.
#   --dataScale X        scale the observed histogram before writing    [1]
#   --noBinByBinStat     drop the template's own statistical nuisance
#   --postfixSuffix S    appended to the output postfix                  [""]
#   --mzDir DIR          base output dir   [/ceph/submit/data/user/p/pmlugato/mz]
#   --tensorDir DIR      [<mzDir>/tensors]     --fitDir DIR  [<mzDir>/fitresults]
#   --logDir DIR         [<mzDir>/logs]        --plotDir DIR [<mzDir>/plots/cvhnano_<postfix>]
#   --image PATH         singularity image
#   --fitExtra "ARGS"    extra arguments appended to the rabbit_fit.py call, e.g.
#                        "--minimizerMethod trust-exact" when trust-krylov stalls
#                        on the observation. Label the run with --postfixSuffix.
#   --dryRun             write the commands file and print the plan; run nothing
#
# Examples:
#   btojpsik_cvhnano_fit.sh $HIST full
#   btojpsik_cvhnano_fit.sh $HIST full --bkgModel bernstein
#   btojpsik_cvhnano_fit.sh $HIST inclusive --observed
#   btojpsik_cvhnano_fit.sh $HIST full -- --pseudoData mcTotal
#   btojpsik_cvhnano_fit.sh $HIST full --bkgModel exppoly --etaRebin 2 \
#       --normAxes fit_bachelor_eta --scaleVar --mcBkg peak \
#       --resVar 0.005 --resAxes fit_bachelor_eta --signalScale 22.1 \
#       --massLow 5.05 --massHigh 5.45 --massBins 20

set -uo pipefail

usage() { sed -n '2,/^set -uo/p' "${BASH_SOURCE[0]}" | sed '$d' | sed 's/^# \{0,1\}//'; }

# Captured before any shift, so the commands file records the whole invocation
# rather than whatever was left after the options were consumed.
INVOCATION="$*"
HIST=${1:?usage: $0 <histmaker.hdf5> <config> [options] [-- extra tensor args...]}
CONFIG=${2:?config must be one of: full ptcharge inclusive}
shift 2

# --- refuse the retired environment knobs ------------------------------------
# Silently ignoring them would be the exact failure mode the flags exist to
# remove, so this is an error, not a warning.
declare -A RETIRED=(
  [BKGMODEL]=--bkgModel [POLYORDER]=--polyOrder [ASIMOV]=--observed
  [SCALEVAR]=--scaleVar [SCALEVARSCALE]=--scaleVarScale
  [SKIPSCALEPARAMS]=--skipScaleParams [SCAN2D]=--scan2D
  [ETA_REBIN]=--etaRebin [PT_BINS]=--ptBins [SLOPE_AXES]=--slopeAxes
  [CURV_AXES]=--curvAxes [NORM_AXES]=--normAxes [MCBKG]=--mcBkg
  [MCBKG_NORM_AXES]=--mcBkgNormAxes [RESVAR]=--resVar [RESCONSTR]=--resConstrained
  [DATASCALE]=--dataScale [NOBBB]=--noBinByBinStat [POSTFIX_SUFFIX]=--postfixSuffix
  [TENSOR_DIR]=--tensorDir [FIT_DIR]=--fitDir [LOG_DIR]=--logDir [PLOTDIR]=--plotDir
)
STALE=()
for name in "${!RETIRED[@]}"; do
  [[ -n ${!name+x} ]] && STALE+=("$name -> ${RETIRED[$name]}")
done
if [[ ${#STALE[@]} -gt 0 ]]; then
  echo "ERROR: environment variables are no longer read by $(basename "$0")." >&2
  echo "The following are set in your shell; unset them and pass the flag instead:" >&2
  printf '  %s\n' "${STALE[@]}" >&2
  exit 2
fi

# --- defaults ------------------------------------------------------------------
WREM_BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
MZ=/ceph/submit/data/user/p/pmlugato/mz
TENSOR_DIR=""; FIT_DIR=""; LOG_DIR=""; PLOTDIR=""
IMAGE=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest

BKGMODEL=exp
# Asimov by default. A failure present on both Asimov and the observation is a
# model or conditioning problem; one present only on the observation is not, and
# that distinction is the whole reason both are run.
ASIMOV=1
ETA_REBIN=1
PT_BINS=3
# Axes the exponential slope varies over. A subset of the cell axes: the
# background slope is driven by acceptance and combinatorial kinematics, both
# primarily eta-dependent, so one slope per eta bin is a physics statement as
# well as a conditioning one.
SLOPE_AXES=fit_bachelor_eta
# Unset means "follow --slopeAxes"; --curvAxes "" means one global curvature,
# which is the measured starting point: the postfit residual is the same shape
# in every cell, so a per-cell curvature would fit 16 copies of one number.
CURV_AXES_SET=0; CURV_AXES=""
POLYORDER=2
MCBKG=none
DATASCALE=1
# Resolution nuisance: smearing width in GeV. ~0.003 from the matched-purity,
# matched-window fit (openspec aem 17.5). The 0.0089 once recommended here was
# sqrt(17.29^2 - 14.81^2), a data-versus-MC-SIGNAL comparison, and 0.005 was
# the wide-window value where the background model absorbs shape.
RESVAR=0
RESCONSTR=0
RES_AXES=""
SIGNAL_SCALE=""
REF_YIELD=""
MASS_LOW=""; MASS_HIGH=""; MASS_BINS=""
BKGSEED=flat
MCBKG_NORM_AXES=""
# Axes the signal normalisation varies over. Empty means every cell axis.
# Setting it to a subset is the first lever to pull if the fit fails on
# statistics, and it is preferable to coarsening eta, which discards the
# resolution structure the channel is measuring.
NORM_AXES=""
# Off by default so the pre-existing no-systematics configuration is still
# exactly reproducible; the campaign turns it on explicitly.
SCALEVAR=0
SCALEVARSCALE=1.0
SKIPSCALEPARAMS=""
# The A-e degeneracy is the one thing the per-nuisance sensitivities cannot see,
# so a 2D scan on a chosen pair is a first-class output rather than a follow-up.
SCAN2D=""
NOBBB=0
POSTFIX_SUFFIX=""
FIT_EXTRA=""
DRYRUN=0

# --- options -------------------------------------------------------------------
need() { [[ $# -ge 2 ]] || { echo "$1 needs a value" >&2; exit 2; }; }
EXTRA=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bkgModel)        need "$@"; BKGMODEL=$2; shift 2 ;;
    --polyOrder)       need "$@"; POLYORDER=$2; shift 2 ;;
    --observed)        ASIMOV=0; shift ;;
    --asimov)          ASIMOV=1; shift ;;
    --scaleVar)        SCALEVAR=1; shift ;;
    --scaleVarScale)   need "$@"; SCALEVARSCALE=$2; shift 2 ;;
    --skipScaleParams) need "$@"; SKIPSCALEPARAMS=$2; shift 2 ;;
    --scan2D)          need "$@"; SCAN2D=$2; shift 2 ;;
    --etaRebin)        need "$@"; ETA_REBIN=$2; shift 2 ;;
    --ptBins)          need "$@"; PT_BINS=$2; shift 2 ;;
    --slopeAxes)       need "$@"; SLOPE_AXES=$2; shift 2 ;;
    --curvAxes)        [[ $# -ge 2 ]] || { echo "--curvAxes needs a value (\"\" for one global curvature)" >&2; exit 2; }
                       CURV_AXES=$2; CURV_AXES_SET=1; shift 2 ;;
    --normAxes)        need "$@"; NORM_AXES=$2; shift 2 ;;
    --mcBkg)           need "$@"; MCBKG=$2; shift 2 ;;
    --mcBkgNormAxes)   need "$@"; MCBKG_NORM_AXES=$2; shift 2 ;;
    --resVar)          need "$@"; RESVAR=$2; shift 2 ;;
    --resConstrained)  RESCONSTR=1; shift ;;
    --resAxes)         need "$@"; RES_AXES=$2; shift 2 ;;
    --signalScale)     need "$@"; SIGNAL_SCALE=$2; shift 2 ;;
    --refYield)        need "$@"; REF_YIELD=$2; shift 2 ;;
    --massLow)         need "$@"; MASS_LOW=$2; shift 2 ;;
    --massHigh)        need "$@"; MASS_HIGH=$2; shift 2 ;;
    --massBins)        need "$@"; MASS_BINS=$2; shift 2 ;;
    --bkgSeed)         need "$@"; BKGSEED=$2; shift 2 ;;
    --dataScale)       need "$@"; DATASCALE=$2; shift 2 ;;
    --noBinByBinStat)  NOBBB=1; shift ;;
    --postfixSuffix)   need "$@"; POSTFIX_SUFFIX=$2; shift 2 ;;
    --mzDir)           need "$@"; MZ=$2; shift 2 ;;
    --tensorDir)       need "$@"; TENSOR_DIR=$2; shift 2 ;;
    --fitDir)          need "$@"; FIT_DIR=$2; shift 2 ;;
    --logDir)          need "$@"; LOG_DIR=$2; shift 2 ;;
    --plotDir)         need "$@"; PLOTDIR=$2; shift 2 ;;
    --image)           need "$@"; IMAGE=$2; shift 2 ;;
    --fitExtra)        need "$@"; FIT_EXTRA=$2; shift 2 ;;
    --dryRun)          DRYRUN=1; shift ;;
    -h|--help)         usage; exit 0 ;;
    --)                shift; EXTRA=("$@"); break ;;
    -*)                echo "unknown option '$1' (extra tensor-writer args go after '--')" >&2; exit 2 ;;
    *)                 echo "unexpected positional '$1' (extra tensor-writer args go after '--')" >&2; exit 2 ;;
  esac
done
[[ $CURV_AXES_SET == 1 ]] || CURV_AXES=$SLOPE_AXES
# A per-eta resolution with no reference width is a request for nothing, and
# silently emitting nothing is the failure mode the flags exist to remove.
if [[ -n $RES_AXES && $RESVAR == 0 ]]; then
  echo "--resAxes '$RES_AXES' needs --resVar <width in GeV> (0.005 for the per-eta primary)" >&2
  exit 2
fi
if [[ -n $MASS_LOW || -n $MASS_HIGH || -n $MASS_BINS ]]; then
  if [[ -z $MASS_LOW || -z $MASS_HIGH || -z $MASS_BINS ]]; then
    echo "--massLow, --massHigh and --massBins must be given together" >&2; exit 2
  fi
fi
TENSOR_DIR=${TENSOR_DIR:-$MZ/tensors}
FIT_DIR=${FIT_DIR:-$MZ/fitresults}
LOG_DIR=${LOG_DIR:-$MZ/logs}

CHANNEL=btojpsik_cvhnano
SIGNAL=BuJpsiK
BKG=BuJpsiKBkg

PT=fit_bachelor_pt
ETA=fit_bachelor_eta
CH=fit_bachelor_charge
MASS=jointCvh_mass

# `inclusive` keeps the charge axis rather than collapsing to a single cell.
# Both param models key their parameters on cell axes, so a one-cell tensor
# would leave them with nothing to key on and force a different model just for
# that configuration. Two charge cells give the same total yield -- it is a sum
# over two numbers instead of one -- and check that the charges agree at the same
# time, so nothing is lost by keeping the axis.
case "$CONFIG" in
  full)      COLLAPSE=() ;;
  ptcharge)  COLLAPSE=("$ETA") ;;
  inclusive) COLLAPSE=("$PT" "$ETA") ;;
  *) echo "unknown config '$CONFIG'; expected full, ptcharge or inclusive" >&2; exit 2 ;;
esac

# The cell axes that survive, in the order the tensor stores them.
CELL_AXES=()
for ax in "$PT" "$ETA" "$CH"; do
  keep=1
  for dropped in ${COLLAPSE+"${COLLAPSE[@]}"}; do
    [[ "$ax" == "$dropped" ]] && keep=0
  done
  [[ $keep == 1 ]] && CELL_AXES+=("$ax")
done
CELL_CSV=$(IFS=,; echo "${CELL_AXES[*]}")

POSTFIX="${CONFIG}_${BKGMODEL}"
if [[ $ASIMOV == 1 ]]; then POSTFIX="${POSTFIX}_asimov"; else POSTFIX="${POSTFIX}_data"; fi
[[ $ETA_REBIN != 1 ]] && POSTFIX="${POSTFIX}_eta${ETA_REBIN}"
[[ -n $NORM_AXES ]] && POSTFIX="${POSTFIX}_normgrouped"
[[ $SCALEVAR == 1 ]] && POSTFIX="${POSTFIX}_aem"
[[ $NOBBB == 1 ]] && POSTFIX="${POSTFIX}_nobbb"
[[ $BKGMODEL == exppoly ]] && POSTFIX="${POSTFIX}_n${POLYORDER}"
[[ $MCBKG != none ]] && POSTFIX="${POSTFIX}_mc${MCBKG}"
[[ $DATASCALE != 1 ]] && POSTFIX="${POSTFIX}_ds$(echo "$DATASCALE" | tr -d '.')"
[[ $RESVAR != 0 ]] && POSTFIX="${POSTFIX}_res$(echo "$RESVAR" | tr -d '.')"
[[ $RESCONSTR == 1 ]] && POSTFIX="${POSTFIX}c"
# _reseta / _resetapt: the axis names with the fit_bachelor_ prefix stripped
# and joined, so the postfix says which axes the resolution varies over.
[[ -n $RES_AXES ]] && POSTFIX="${POSTFIX}_res$(echo "$RES_AXES" | sed 's/fit_bachelor_//g; s/,//g')"
[[ -n $SIGNAL_SCALE ]] && POSTFIX="${POSTFIX}_ss$(echo "$SIGNAL_SCALE" | tr -d '.')"
[[ -n $MASS_BINS ]] && POSTFIX="${POSTFIX}_m${MASS_BINS}"
[[ $BKGSEED != flat ]] && POSTFIX="${POSTFIX}_seed${BKGSEED}"
# Dots out of the label: the postfix becomes a filename, and "x10.0" produced a
# tensor the fit step then could not find.
[[ $SCALEVARSCALE != 1.0 ]] && POSTFIX="${POSTFIX}_x${SCALEVARSCALE//./p}"
[[ -n $SKIPSCALEPARAMS ]] && POSTFIX="${POSTFIX}_no$(echo "$SKIPSCALEPARAMS" | tr -d ' ')"
POSTFIX="${POSTFIX}${POSTFIX_SUFFIX}"

TENSOR="$TENSOR_DIR/btojpsik_cvhnano_tensor_${POSTFIX}.hdf5"
FIT="$FIT_DIR/cvhnano_${POSTFIX}_fit.hdf5"
PLOTREADY="$FIT_DIR/cvhnano_${POSTFIX}_plotready.hdf5"
PLOTDIR=${PLOTDIR:-$MZ/plots/cvhnano_${POSTFIX}}

mkdir -p "$TENSOR_DIR" "$FIT_DIR" "$LOG_DIR"

run() {
  singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc \
    "cd '$WREM_BASE' && source setup.sh >/dev/null 2>&1 && $1"
}

# --- tensor -----------------------------------------------------------------
TENSOR_CMD="/opt/venv/bin/python scripts/rabbit/btojpsik_cvhnano_tensor.py \
  -i '$HIST' -o '$TENSOR_DIR' -p '$POSTFIX' \
  --ptBins $PT_BINS --etaRebin $ETA_REBIN"
[[ ${#COLLAPSE[@]} -gt 0 ]] && TENSOR_CMD="$TENSOR_CMD --collapse ${COLLAPSE[*]}"
if [[ $SCALEVAR == 1 ]]; then
  TENSOR_CMD="$TENSOR_CMD --scaleVariations --scaleVarScale $SCALEVARSCALE"
  [[ -n $SKIPSCALEPARAMS ]] && TENSOR_CMD="$TENSOR_CMD --skipScaleParams $SKIPSCALEPARAMS"
fi
[[ $MCBKG != none ]] && TENSOR_CMD="$TENSOR_CMD --mcBkgTemplate $MCBKG"
[[ $DATASCALE != 1 ]] && TENSOR_CMD="$TENSOR_CMD --dataScale $DATASCALE"
[[ $RESVAR != 0 ]] && TENSOR_CMD="$TENSOR_CMD --resolutionVariation $RESVAR"
[[ $RESCONSTR == 1 ]] && TENSOR_CMD="$TENSOR_CMD --resolutionConstrained"
[[ -n $RES_AXES ]] && TENSOR_CMD="$TENSOR_CMD --resolutionAxes ${RES_AXES//,/ }"
[[ -n $SIGNAL_SCALE ]] && TENSOR_CMD="$TENSOR_CMD --signalScale $SIGNAL_SCALE"
[[ -n $REF_YIELD ]] && TENSOR_CMD="$TENSOR_CMD --signalScaleReferenceYield $REF_YIELD"
[[ -n $MASS_BINS ]] && TENSOR_CMD="$TENSOR_CMD --massLow $MASS_LOW --massHigh $MASS_HIGH --massBins $MASS_BINS"
[[ $BKGSEED != flat ]] && TENSOR_CMD="$TENSOR_CMD --bkgSeed $BKGSEED"
[[ ${#EXTRA[@]} -gt 0 ]] && TENSOR_CMD="$TENSOR_CMD ${EXTRA[*]}"

# --- param models -----------------------------------------------------------
# Two composed models. rabbit_fit.py wraps several --paramModel occurrences in
# CompositeParamModel by itself, so no wrapper class is needed.
NORM_CSV=${NORM_AXES:-$CELL_CSV}
MODELS="--paramModel AxisNormModel $CHANNEL $SIGNAL $NORM_CSV"
# The MC template normalisations float too, on their own grouping. They are
# seeded at the simulation's own prediction (a fitted norm of 1 reproduces it),
# so the fitted value is a measurement of the composition, not of a scale.
MCBKG_NORM_CSV=${MCBKG_NORM_AXES:-$NORM_CSV}
case $MCBKG in
  all)   MCPROCS="BuJpsiKMCBkg" ;;
  comb)  MCPROCS="BuJpsiKMCComb" ;;
  split) MCPROCS="BuJpsiKMCComb,BuJpsiKMCPeak" ;;
  peak)  MCPROCS="BuJpsiKMCPeak" ;;
  none)  MCPROCS="" ;;
  *) echo "unknown --mcBkg '$MCBKG'; expected none, all, comb, split or peak" >&2; exit 2 ;;
esac
[[ -n $MCPROCS ]] && MODELS="$MODELS --paramModel AxisNormModel $CHANNEL $MCPROCS $MCBKG_NORM_CSV"
# Keep only slope axes that survived the collapse; naming a collapsed axis
# would be rejected by the model rather than ignored.
SLOPE_CSV=""
for ax in ${SLOPE_AXES//,/ }; do
  case ",$CELL_CSV," in *",$ax,"*) SLOPE_CSV="${SLOPE_CSV:+$SLOPE_CSV,}$ax" ;; esac
done
if [[ $BKGMODEL == exp ]]; then
  MODELS="$MODELS --paramModel AxisExpModel $CHANNEL $BKG $MASS $CELL_CSV${SLOPE_CSV:+ $SLOPE_CSV}"
elif [[ $BKGMODEL == expquad ]]; then
  # Alias for exppoly order 2. AxisExpQuadModel was lost when the rabbit
  # submodule was reset to havyn/upstream-main-havyn on 2026-08-27; order 2 of
  # the poly model was verified equivalent to it (span residual 3.2e-07, and the
  # two Asimov fits agree to 2.2338972927e-29 vs 2.2338976236e-29), so the alias
  # loses nothing and there is one class to maintain instead of two.
  [[ -z $SLOPE_CSV ]] && SLOPE_CSV="$CELL_CSV"
  MODELS="$MODELS --paramModel AxisExpPolyModel $CHANNEL $BKG $MASS $CELL_CSV 2 $SLOPE_CSV"
elif [[ $BKGMODEL == exppoly ]]; then
  # exp(lnAmpl + sum_k c_k P_k), P_k orthogonalised on the actual bin grid by
  # Gram-Schmidt against the constant and each lower order. Positive for any real
  # parameters, and every order nests the one below at an interior point, so
  # adjacent orders differ by a likelihood ratio for one parameter.
  [[ -z $SLOPE_CSV ]] && SLOPE_CSV="$CELL_CSV"
  # One grouping given; the model reuses it for every higher order.
  MODELS="$MODELS --paramModel AxisExpPolyModel $CHANNEL $BKG $MASS $CELL_CSV $POLYORDER $SLOPE_CSV"
elif [[ $BKGMODEL == bernstein ]]; then
  # NB first order only: c0*(1-x) + c1*x. This is a LINEAR background with no
  # curvature -- it is not the curved alternative to exp, expquad is.
  MODELS="$MODELS --paramModel AxisBernsteinModel $CHANNEL $BKG $MASS $CELL_CSV"
else
  echo "unknown --bkgModel '$BKGMODEL'; expected exp, expquad, exppoly or bernstein" >&2; exit 2
fi

# --computeHistCov only on Asimov. On the observation the per-cell projection
# (16 cells x 80 mass bins = 1280 bins) has a rank-deficient prediction
# covariance, and the chi-square that consumes it dies with "Input matrix is not
# invertible" *after* the fit has converged -- taking the results file with it.
# The fit itself is unaffected; only the goodness-of-fit bookkeeping is.
HISTCOV=""
if [[ $ASIMOV == 1 ]]; then
  HISTCOV="--computeHistCov"
else
  # --noChi2 as well, and this one is not optional. On the observation the
  # goodness-of-fit routine inverts a prediction covariance that is
  # rank-deficient over the 1280-bin per-cell projection, dies with "Input
  # matrix is not invertible" *after* the fit has converged, and takes the
  # results file with it -- the output is left holding only 'meta'. The fit is
  # unaffected; only the chi-square bookkeeping is, and it has to be computed
  # separately.
  HISTCOV="--noChi2"
fi

FIT_CMD="/opt/venv/bin/python rabbit/bin/rabbit_fit.py '$TENSOR' -v 4 --diagnostics \
  -o '$FIT_DIR' --outname 'cvhnano_${POSTFIX}_fit.hdf5' \
  --saveHists --saveHistsPerProcess --computeHistErrors \
  --computeHistErrorsPerProcess $HISTCOV \
  $MODELS \
  -m Project $CHANNEL $MASS \
  -m Project $CHANNEL ${CELL_AXES[*]} $MASS"
if [[ $NOBBB == 1 ]]; then
  # Turning off the template's own statistical nuisance is a diagnostic, not a
  # configuration: it says how much of the constraint the simulated sample size
  # is costing, which is a different question from how much the data costs.
  FIT_CMD="$FIT_CMD --noBinByBinStat"
fi
if [[ $SCALEVAR == 1 ]]; then
  # Impacts and the nuisance covariance are the point of adding the nuisances:
  # without --doImpacts the fit reports uncertainties and never says what they
  # are made of, and without --computeHistCov the A-e correlation cannot be read
  # off at all.
  FIT_CMD="$FIT_CMD --doImpacts --globalImpacts"
  # --chisqFit only on Asimov. On the observation it raises
  #   "Bins in 'nobs <= 0' encountered, chi^2 fit can not be performed"
  # because real data has empty bins in the sparse cells, and the run dies after
  # writing a results-less output file. The Poisson likelihood is the default and
  # the right choice there.
  [[ $ASIMOV == 1 ]] && FIT_CMD="$FIT_CMD --chisqFit"
fi
if [[ -n $SCAN2D ]]; then
  FIT_CMD="$FIT_CMD --scan2D $SCAN2D --scanPoints 21 --scanRange 3"
fi
# rabbit's -t is the opposite way round from combinetf's: -1 is the Asimov toy
# (and the default), 0 fits the observation. Getting this backwards mislabels
# every result in the campaign, and the "does it fail on Asimov too" attribution
# depends on it entirely.
# --unblind on the observation, for the fit and (below) the plot-ready pass:
# rabbit's default blinding multiplies every POI by exp(N(0, 5)), which for a
# calibration channel is a random start up to e^10 off in every cell, and the
# externalPostfit pass then evaluates the model at blinded parameters and dies
# on a NaN Hessian.
if [[ $ASIMOV == 1 ]]; then FIT_CMD="$FIT_CMD -t -1"; else FIT_CMD="$FIT_CMD -t 0 --unblind"; fi
[[ -n $FIT_EXTRA ]] && FIT_CMD="$FIT_CMD $FIT_EXTRA"

PLOT_CMD="/opt/venv/bin/python rabbit/bin/rabbit_fit.py '$TENSOR' \
  --externalPostfit '$FIT' --saveHists --computeHistErrors \
  -o '$FIT_DIR' --outname 'cvhnano_${POSTFIX}_plotready.hdf5' \
  $MODELS \
  -m Project $CHANNEL $MASS \
  -m Project $CHANNEL ${CELL_AXES[*]} $MASS"
if [[ $ASIMOV == 1 ]]; then PLOT_CMD="$PLOT_CMD -t -1"; else PLOT_CMD="$PLOT_CMD -t 0 --unblind"; fi

COMMANDS="$LOG_DIR/cvhnano_${POSTFIX}_commands.sh"
{
  echo "#!/usr/bin/env bash"
  echo "# written by $(basename "$0"); run these inside the container"
  echo "# invocation: $(basename "$0") $INVOCATION"
  echo "$TENSOR_CMD"
  echo "$FIT_CMD"
  echo "$PLOT_CMD"
} > "$COMMANDS"
chmod +x "$COMMANDS"

echo "=== $CONFIG / $BKGMODEL / $([[ $ASIMOV == 1 ]] && echo Asimov || echo observed) ==="
echo "cell axes: ${CELL_AXES[*]}"
echo "tensor:    $TENSOR"
echo "fit:       $FIT"
echo "commands:  $COMMANDS"
echo
if [[ $DRYRUN == 1 ]]; then
  echo "--dryRun: nothing executed. Commands:"
  cat "$COMMANDS"
  exit 0
fi

run "$TENSOR_CMD" 2>&1 | tee "$LOG_DIR/cvhnano_${POSTFIX}_tensor.log"
if [[ ${PIPESTATUS[0]} != 0 ]]; then
  echo "TENSOR FAILED; see $LOG_DIR/cvhnano_${POSTFIX}_tensor.log" >&2
  exit 1
fi

run "$FIT_CMD" 2>&1 | tee "$LOG_DIR/cvhnano_${POSTFIX}_fit.log"
if [[ ${PIPESTATUS[0]} != 0 ]]; then
  # Not swallowed, and not a reason to stop the campaign: recorded, with the log
  # left for the attribution step, which needs the failures as much as the
  # successes.
  echo >&2
  echo "FIT DID NOT COMPLETE. Log: $LOG_DIR/cvhnano_${POSTFIX}_fit.log" >&2
  echo "Attribute it before reporting: which parameter sat at a boundary, which" >&2
  echo "cell was empty, which pair was degenerate, or whether it is statistics --" >&2
  echo "and whether the same failure appears on Asimov." >&2
  exit 0
fi

# --- convergence gate -------------------------------------------------------
# A zero exit status is not convergence. rabbit writes its results whenever the
# Hessian passes the Cholesky factorisation, and a Hessian at a NON-STATIONARY
# point factorises perfectly well while being no covariance at all. The
# normgrouped exp data fit ran to iteration 7199 with edmval 72668, wrote a
# 442 kB results file with plots and impacts, and every data-side number taken
# from it -- postfit chi-square, correlations, pulls, impacts -- was meaningless.
#
# **edmval is the only verdict.** Two things that look like verdicts are not:
#
#   scipy's "success" flag: trust-krylov reports success: False on fits that have
#   converged to edmval 1e-17. Measured on three toys, all three False. Using it
#   as a verdict rejects every good fit, so it is a warning here and nothing more.
#
#   "Maximum number of iterations has been exceeded": means the minimiser ran out
#   of budget, which usually shows up in edmval anyway. Surfaced loudly, because a
#   budget-limited fit is worth knowing about, but not itself fatal.
#
# Negative edmval fails. A failed fit elsewhere in this analysis reported
# edmvalref = -0.0, and "-0.0 <= 1e-5" is true, so the test is two-sided.
# Multiple edmvals (a multi-toy run prints one per toy) are all checked, and the
# worst decides.
FITLOG="$LOG_DIR/cvhnano_${POSTFIX}_fit.log"
EDMVALS=$(grep -oE 'edmval: [0-9.eE+-]+' "$FITLOG" | awk '{print $2}')
if grep -q 'Maximum number of iterations has been exceeded' "$FITLOG"; then
  echo "WARNING: minimiser hit its iteration cap; consider --minimizerMaxiter." >&2
fi
if grep -q 'success: False' "$FITLOG"; then
  echo "note: scipy reported success: False, which trust-krylov also does on" >&2
  echo "converged fits. edmval decides, not this." >&2
fi
if [[ -z $EDMVALS ]]; then
  echo "*** $POSTFIX: no edmval in $FITLOG -- convergence UNVERIFIED ***" >&2
  echo "Outputs are diagnostics, not results." >&2
else
  BAD=$(echo "$EDMVALS" | awk 'BEGIN{n=0} {if (!($1+0 <= 1e-5 && $1+0 >= 0)) n++} END{print n}')
  WORST=$(echo "$EDMVALS" | awk 'BEGIN{w=0} {v=$1+0; if (v<0) v=-v; if (v>w) w=v} END{printf "%g", w}')
  NTOT=$(echo "$EDMVALS" | wc -l)
  if [[ $BAD -gt 0 ]]; then
    echo "*** $POSTFIX DID NOT CONVERGE: $BAD of $NTOT edmval(s) outside" >&2
    echo "[0, 1e-5], worst |edmval| = $WORST. A Hessian evaluated away from the" >&2
    echo "minimum is not a covariance: correlations, pulls and impacts from this" >&2
    echo "fit are not results. Outputs below are diagnostics only. ***" >&2
  else
    echo "converged: $NTOT edmval(s), worst |edmval| = $WORST"
  fi
fi

run "$PLOT_CMD" 2>&1 | tee "$LOG_DIR/cvhnano_${POSTFIX}_plotready.log"

# --- readout ----------------------------------------------------------------
# rabbit's own tools, deliberately, rather than bespoke plotting: they read the
# fit file's own conventions for prefit/postfit, impacts and the nuisance
# covariance, so there is nothing to get wrong in re-deriving them.
mkdir -p "$PLOTDIR"

readout() {
  local label=$1; shift
  echo "--- $label"
  run "$1" 2>&1 | tee "$LOG_DIR/cvhnano_${POSTFIX}_${label}.log" | tail -5
  if [[ ${PIPESTATUS[0]} != 0 ]]; then
    echo "   $label did not complete; see $LOG_DIR/cvhnano_${POSTFIX}_${label}.log" >&2
  fi
}

# The mapping has to be named explicitly. Without -m the script loops over "all
# by default", recognises nothing, and exits 0 having written no plot at all --
# a silent no-op, which is why it is spelled out here.
MAP_MASS="-m Project $CHANNEL $MASS"
MAP_CELLS="-m Project $CHANNEL ${CELL_AXES[*]} $MASS"

readout hists_prefit "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists.py '$PLOTREADY' -o '$PLOTDIR' --prefit --dataName Data $MAP_MASS"
readout hists_postfit "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists.py '$PLOTREADY' -o '$PLOTDIR' --dataName Data $MAP_MASS"
readout hists_cells "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists.py '$PLOTREADY' -o '$PLOTDIR' --dataName Data $MAP_CELLS --selectionAxes ${CELL_AXES[*]}"
readout uncertainties "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists_uncertainties.py '$PLOTREADY' -o '$PLOTDIR' $MAP_MASS"

if [[ $SCALEVAR == 1 ]]; then
  readout impacts "/opt/venv/bin/python rabbit/bin/rabbit_plot_pulls_and_impacts.py '$FIT' -o '$PLOTDIR' --mode both --showNumbers --otherExtensions png"
  readout impacts_text "/opt/venv/bin/python rabbit/bin/rabbit_print_impacts.py '$FIT'"
  # --correlation, not the raw covariance: the question is whether A and e are
  # degenerate, and a covariance in mixed units cannot answer it. --params has to
  # name the nuisances, or the matrix is dominated by the normalisations and the
  # block that matters is invisible.
  # The eta-group count is read from the tensor log rather than assumed: this
  # once hard-coded groups 0-3 and silently plotted 12 of 72 nuisances.
  TENSORLOG="$LOG_DIR/cvhnano_${POSTFIX}_tensor.log"
  NGROUPS=$(grep -oE '^.*scaleA[0-9]+ ' "$TENSORLOG" | grep -oE 'scaleA[0-9]+' | sort -u | wc -l)
  [[ $NGROUPS -gt 0 ]] || NGROUPS=4
  SCALE_PARAMS=$(for prm in A e M; do
    case " $SKIPSCALEPARAMS " in *" $prm "*) continue ;; esac
    for ((g=0; g<NGROUPS; g++)); do printf 'scale%s%s ' "$prm" "$g"; done
  done)
  # No --showNumbers: the printed values overlap into an unreadable smear, and
  # the question the plot answers -- is the matrix diagonal -- is answered by
  # the colours.
  readout correlation "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -o '$PLOTDIR' --correlation --params $SCALE_PARAMS"
  # Scale against resolution: the three scale nuisances of the middle eta group
  # together with EVERY resolution nuisance, so the plot shows which resolution
  # bin talks to which scale group without assuming the two eta indexings agree.
  RES_PARAMS=$(grep -oE 'sigRes[A-Za-z0-9_]*' "$TENSORLOG" | sort -u | tr '\n' ' ')
  if [[ -n $RES_PARAMS ]]; then
    MID=$((NGROUPS / 2))
    MID_SCALE=$(for prm in A e M; do
      case " $SKIPSCALEPARAMS " in *" $prm "*) continue ;; esac
      printf 'scale%s%s ' "$prm" "$MID"
    done)
    readout correlation_res "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -o '$PLOTDIR' --correlation --postfix scaleVsRes --params $MID_SCALE $RES_PARAMS"
  fi
fi
if [[ -n $SCAN2D ]]; then
  readout scan2d "/opt/venv/bin/python rabbit/bin/rabbit_plot_likelihood_scan2D.py '$FIT'     -o '$PLOTDIR'"
fi

echo
echo "done: $FIT"
echo "      $PLOTREADY"
echo "      $PLOTDIR"
