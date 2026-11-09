#!/usr/bin/env bash
# Post-processing of the B+ statistics scaling test (notes/btojpsik_stat_scaling_test.md):
# overlay impact plots, the summary table, and an index page. Idempotent; re-run
# after any chain finishes.
set -uo pipefail
WREM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FITS=/ceph/submit/data/user/p/pmlugato/mz/fitresults
OUT=$HOME/public_html/mz/joint_stat_scaling_260916
IMAGE=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest
run() { singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc "cd '$WREM' && source setup.sh >/dev/null 2>&1 && $1" 2>&1 | grep -v 'absl\|TensorFlow\|cpu_feature\|rebuild\|FutureWarning\|metadata was not set'; }
fit() { echo "$FITS/joint_mz_jpsi_bplus_$1_fit.hdf5"; }
mkdir -p "$OUT/compare"

COMMON="--poi massShiftZ100MeV --mode group --oneSidedImpacts --showNumbers --noPulls --scaleImpacts 100 --impactTitle 'Impact on m<sub>Z</sub> [MeV]' --config scripts/rabbit/mz_impacts_plotconfig.py -o $OUT/compare --otherExtensions png pdf"
overlay() { # <fit tag> <ref tag> <name> <refName> <grouping> <postfix>
  [[ -f $(fit $1) && -f $(fit $2) ]] || { echo "skip $6 (missing fit)"; return; }
  run "python rabbit/bin/rabbit_plot_pulls_and_impacts.py $(fit $1) -r $(fit $2) --name '$3' --refName '$4' --grouping $5 --postfix $6 $COMMON" | grep -i "write output file.*png\|error\|traceback"
}
# full fit: scale groups are named A/e/M in every full fit -> mzscale_aem; Z + B+: scaleA/e/M -> mzscale_scale
for N in 1 10 50; do
  T=statscale_mc2x_data$N
  overlay $T res5_260914            "MC 2x data, data x$N"  "MC as produced, data x1 (res5)" mzscale_aem   full_mc2x_data${N}_vs_res5
  overlay ${T}_noJpsi res5_260914_noJpsi "Z+B: MC 2x data, data x$N" "Z+B: MC as produced (res5)" mzscale_scale ZB_mc2x_data${N}_vs_res5
  [[ $N != 1 ]] && overlay $T statscale_mc2x_data1 "data x$N" "data x1 (nominal)" mzscale_aem full_mc2x_data${N}_vs_nominal
  [[ $N != 1 ]] && overlay ${T}_noJpsi statscale_mc2x_data1_noJpsi "Z+B: data x$N" "Z+B: data x1 (nominal)" mzscale_scale ZB_mc2x_data${N}_vs_nominal
done
# full vs Z+B at each point
for N in 1 10 50; do
  # the Z + B+ fit only has the scaleA/e/M group names; the full fit has both, so use the scale naming here
  overlay statscale_mc2x_data$N statscale_mc2x_data${N}_noJpsi "Z + J/ψ + B⁺" "Z + B⁺" mzscale_scale full_vs_ZB_mc2x_data$N
done

SPECS=""
for spec in "res5 (MC as produced)=res5_260914" "res5 Z+B=res5_260914_noJpsi" "Z+J/psi only=res5_260914_noB" \
            "MC2x data1=statscale_mc2x_data1" "MC2x data1 Z+B=statscale_mc2x_data1_noJpsi" \
            "MC2x data10=statscale_mc2x_data10" "MC2x data10 Z+B=statscale_mc2x_data10_noJpsi" \
            "MC2x data50=statscale_mc2x_data50" "MC2x data50 Z+B=statscale_mc2x_data50_noJpsi"; do
  SPECS="$SPECS '${spec%%=*}=$(fit ${spec##*=})'"
done
run "python scripts/rabbit/statscale_report.py --fits $SPECS -o $OUT/statscale_report.md --plot $OUT/statscale_summary.png" | tail -30

python3 - "$OUT" <<'PY'
import os, sys, glob, html
out=sys.argv[1]
L=["<html><head><meta charset='utf-8'><title>B+ statistics scaling test</title><style>body{font-family:sans-serif;max-width:1500px;margin:auto} img{max-width:100%;border:1px solid #ddd} figcaption{color:#52514e;font-size:14px} pre{background:#f6f6f4;padding:12px;white-space:pre-wrap}</style></head><body>",
   "<h1>B+ statistics scaling test (2026-09-16)</h1>",
   "<p>Explanation and log: <code>notes/btojpsik_stat_scaling_test.md</code> in WRemnants. Per-fit rabbit plots (prefit/postfit, pulls, correlations) are in the per-tag subdirectories below.</p>"]
if os.path.exists(f"{out}/statscale_report.md"):
    L.append("<h2>Summary</h2><pre>"+html.escape(open(f"{out}/statscale_report.md").read())+"</pre>")
    if os.path.exists(f"{out}/statscale_summary.png"): L.append("<figure><img src='statscale_summary.png'></figure>")
L.append("<h2>Overlay impact plots</h2>")
for f in sorted(glob.glob(f"{out}/compare/*.png")):
    b=os.path.basename(f); L.append(f"<figure><img src='compare/{b}'><figcaption>{b}</figcaption></figure>")
L.append("<h2>Per-fit directories</h2><ul>")
for d in sorted(p for p in glob.glob(f"{out}/*") if os.path.isdir(p) and os.path.basename(p)!="compare"):
    b=os.path.basename(d); L.append(f"<li><a href='{b}/'>{b}</a> ({len(glob.glob(d+'/*.png'))} png)</li>")
L.append("</ul></body></html>")
open(f"{out}/index.html","w").write("\n".join(L)); print("index written")
PY
