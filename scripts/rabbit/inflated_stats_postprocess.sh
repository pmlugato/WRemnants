#!/usr/bin/env bash
# J/psi + B+ with mB floating, matched unit, Asimov: B+ statistics inflated.
#   as produced : MC as produced, data x1         (joint_mz_jpsi_bplus_matched_mB_noZ_jpsimu_jpsifull)
#   mc2x data1  : MC 2x data, data x1             (--dataScale 1  --mcStatScale 36.4866)
#   mc2x data10 : MC 2x data, data x10            (--dataScale 10 --mcStatScale 364.866)
# Only the B+ statistics change; the J/psi tensor is the same in all three.
# Writes ~/public_html/mz/combined_fits_and_mB_prefit_scale_uncs_matched/07_inflated_stats/.
set -uo pipefail
WREM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$WREM"
F=/ceph/submit/data/user/p/pmlugato/mz/fitresults
P=$HOME/public_html/mz/combined_fits_and_mB_prefit_scale_uncs_matched/07_inflated_stats
IMAGE=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest
run() { singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc "cd '$WREM' && source setup.sh >/dev/null 2>&1 && export OPENBLAS_NUM_THREADS=1 && $1" 2>&1 | grep -v "absl\|TensorFlow\|cpu_feature\|rebuild\|FutureWarning\|metadata was not\|alpha stage\|warnings.warn\|cuInit"; }
have() { [[ -r "$1" ]] || { echo "  skip, missing $(basename $1)"; return 1; }; }
mkdir -p $P/{impacts,constraints,correlations}

REF=$F/joint_mz_jpsi_bplus_matched_mB_noZ_jpsimu_jpsifull_fit.hdf5
D1=$F/joint_mz_jpsi_bplus_matched_mB_mc2x_data1_noZ_jpsimu_jpsifull_fit.hdf5
D10=$F/joint_mz_jpsi_bplus_matched_mB_mc2x_data10_noZ_jpsimu_jpsifull_fit.hdf5
JPA=$F/realmc_D_jpsionly_asimov_fit.hdf5
OLD=$F/realmc_C_asimov_fit.hdf5
for f in $REF $D1 $D10; do have $f || exit 1; done

echo "=== numbers: sigma(mB), conditional chain, group impacts, sigma(A/e/M), rho"
run "python scripts/rabbit/statscale_report.py --poi massShiftB5MeV --poiScale 5 --poiLabel mB \
      --fits 'old unit, as produced=$OLD' 'matched, as produced=$REF' 'matched, MC 2x, data x1=$D1' \
             'matched, MC 2x, data x10=$D10' -o $P/mB_report.md"
run "python scripts/rabbit/matched_units_summary_plot.py -o $P/mB_summary.png \
      --poi massShiftB5MeV --poiScale 5 --poiLabel mB \
      --title 'J/psi + B+ Asimov, mB floating, matched unit: B+ statistics inflated' \
      --fit ref asimov 'J/psi alone' $JPA --fit old asimov 'old unit, as produced' $OLD \
      --fit matched asimov 'matched, as produced' $REF --fit matched asimov 'matched, MC 2x, data x1' $D1 \
      --fit matched asimov 'matched, MC 2x, data x10' $D10"

CFG="--config scripts/rabbit/mz_impacts_plotconfig.py --poi massShiftB5MeV --scaleImpacts 5 --impactTitle 'Impact on m<sub>B</sub> [MeV]' --otherExtensions png pdf"
PI="python rabbit/bin/rabbit_plot_pulls_and_impacts.py"
echo "=== one-sided impacts (grouped, and top 40), each against the as-produced fit"
for s in "data1:$D1:MC 2x, data x1" "data10:$D10:MC 2x, data x10"; do
  IFS=: read -r tag f lab <<< "$s"
  run "$PI $f -r $REF --name '$lab' --refName 'as produced' --mode group --grouping mb_aem \
        --oneSidedImpacts --showNumbers --noPulls --postfix grouped_${tag}_vs_asproduced -o $P/impacts $CFG | grep -c png"
  run "$PI $f -r $REF --name '$lab' --refName 'as produced' --mode ungrouped -n 40 \
        --oneSidedImpacts --postfix top40_${tag}_vs_asproduced -o $P/impacts $CFG | grep -c png"
done
run "$PI $D10 -r $D1 --name 'MC 2x, data x10' --refName 'MC 2x, data x1' --mode group --grouping mb_aem \
      --oneSidedImpacts --showNumbers --noPulls --postfix grouped_data10_vs_data1 -o $P/impacts $CFG | grep -c png"

echo "=== constraints: A_i, e_i, M_i (Asimov, so pulls are zero and the box is the constraint)"
for p in A e M; do
  for s in "data1:$D1:MC 2x, data x1" "data10:$D10:MC 2x, data x10"; do
    IFS=: read -r tag f lab <<< "$s"
    run "$PI $f -r $REF --name '$lab' --refName 'as produced' --mode ungrouped --filters '^${p}[0-9]+\$' \
          -s label -n 24 --oneSidedImpacts --postfix ${p}_${tag}_vs_asproduced -o $P/constraints $CFG | grep -c png"
  done
done
run "python scripts/rabbit/aem_sigma_per_group.py --poi massShiftB5MeV --poiScale 5 --poiLabel mB \
      --fits 'J/psi alone=$JPA' 'matched, as produced=$REF' 'matched, MC 2x, data x1=$D1' 'matched, MC 2x, data x10=$D10' \
      -o $P/constraints/aem_sigma_per_group.png --table $P/constraints/table_aem.md \
      --title 'J/psi + B+ Asimov: per-eta sigma(A/e/M), B+ statistics inflated'"

echo "=== correlations"
for s in "asproduced:$REF" "data1:$D1" "data10:$D10"; do
  tag=${s%%:*}; f=${s#*:}
  for reg in 'barrel_eta12 12' 'endcap_eta0 0'; do set -- $reg
    run "python rabbit/bin/rabbit_plot_cov.py $f -o $P/correlations --correlation --showNumbers --title '' \
          --config scripts/rabbit/mz_impacts_plotconfig.py --postfix corr4x4_${tag}_$1 --params massShiftB5MeV A$2 e$2 M$2 >/dev/null"
  done
  run "python scripts/rabbit/poi_correlation_families.py $f -o $P/correlations/corr_mB_families_${tag}.png \
        --table $P/correlations/corr_mB_families_${tag}.md"
done
echo "=== done"
