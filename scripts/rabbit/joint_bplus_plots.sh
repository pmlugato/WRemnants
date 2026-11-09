#!/usr/bin/env bash
# Plot stage of joint_with_btojpsik.sh, as a separate script so it can be re-run
# on an existing fit without touching the driver while a chain is running.
#
# Usage: joint_bplus_plots.sh <fit.hdf5> <outdir> <nEta> <tag> [--noJpsi] [--noB]
#
# rabbit's own tools throughout. One rabbit_plot_hists call per channel with
# --filterProcs, because the merged tensor's process list is global and a
# single call would legend every channel's processes on every panel (eight
# entries on the Z plot, covering the peak's flank). --titlePos 0 keeps the
# title off the panel; the legend sits upper left, one column, at most four
# entries, clear of every peak.
set -uo pipefail
FIT=$1; PLOTS=$2; N_ETA=$3; TAG=$4; shift 4
NOJPSI=0; NOB=0; NOZ=0; POI=massShiftZ100MeV; POISCALE=100; DATANAME=Asimov
while [[ $# -gt 0 ]]; do case $1 in --noJpsi) NOJPSI=1; shift ;; --noB) NOB=1; shift ;; --noZ) NOZ=1; shift ;; --data) DATANAME=Data; shift ;; --poi) POI=$2; shift 2 ;; --poiScale) POISCALE=$2; shift 2 ;; *) echo "unknown option $1" >&2; exit 2 ;; esac; done
WREM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGE=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest
LOG=/ceph/submit/data/user/p/pmlugato/mz/logs/joint_bplus_plots_${TAG}.log
run() { singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc "cd '$WREM' && source setup.sh >/dev/null 2>&1 && $1"; }
mkdir -p "$PLOTS"; : > "$LOG"
STYLE="--titlePos 0 --legPos 'upper left' --legCols 1 --lowerLegPos 'upper left' --lowerLegCols 1 --dataName $DATANAME"

declare -a CH=(); [[ $NOZ == 0 ]] && CH+=("Z mll Zmumu,Ztautau,Other,PhotonInduced")
[[ $NOJPSI == 0 ]] && CH+=("JPsi_prompt mass sig_JPsi_prompt,background_JPsi_prompt" "JPsi_displaced mass sig_JPsi_displaced,background_JPsi_displaced")
[[ $NOB == 0 ]] && CH+=("btojpsik_cvhnano jointCvh_mass BuJpsiK,BuJpsiKMCPeak,BuJpsiKBkg")
for spec in "${CH[@]}"; do
  set -- $spec; ch=$1; ax=$2; procs=${3//,/ }
  for FT in "--prefit" ""; do
    run "/opt/venv/bin/python rabbit/bin/rabbit_plot_hists.py '$FIT' -o '$PLOTS' -m Project $ch $ax $FT $STYLE --filterProcs $procs" >> "$LOG" 2>&1
  done
done
if [[ -n $POI ]]; then
run "/opt/venv/bin/python rabbit/bin/rabbit_plot_pulls_and_impacts.py '$FIT' -o '$PLOTS' --poi $POI --mode both --num 25 --otherExtensions png --postfix $POI" >> "$LOG" 2>&1
fi
run "/opt/venv/bin/python rabbit/bin/rabbit_print_impacts.py '$FIT'" > "${LOG%.log}_impacts.txt" 2>&1
MID=$((N_ETA / 2)); CLO=$((N_ETA / 2 - 4)); CHI=$((N_ETA / 2 + 3))
run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -o '$PLOTS' --correlation --showNumbers --title '' --config scripts/rabbit/mz_impacts_plotconfig.py --postfix barrel_eta$MID --params $POI A$MID e$MID M$MID" >> "$LOG" 2>&1
run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -o '$PLOTS' --correlation --showNumbers --title '' --config scripts/rabbit/mz_impacts_plotconfig.py --postfix endcap_eta0 --params $POI A0 e0 M0" >> "$LOG" 2>&1
ALLSCALE=$(for p in A e M; do for ((i=0; i<N_ETA; i++)); do printf '%s%s ' "$p" "$i"; done; done)
run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -o '$PLOTS' --correlation --postfix scale_all --params $ALLSCALE" >> "$LOG" 2>&1
CENTRAL=$(for i in $(seq $CLO $CHI); do printf 'A%s e%s M%s ' "$i" "$i" "$i"; done)
run "/opt/venv/bin/python rabbit/bin/rabbit_plot_cov.py '$FIT' -o '$PLOTS' --correlation --postfix poi_central --params $POI $CENTRAL" >> "$LOG" 2>&1
run "/opt/venv/bin/python scripts/plotting/joint_mz_jpsi_summary.py '$FIT' -o '$PLOTS' --tag $TAG" >> "$LOG" 2>&1
echo "plots: $PLOTS ($(ls "$PLOTS" | grep -c '\.png$') png); errors: $(grep -ci 'traceback' "$LOG")"
