#!/usr/bin/env bash
# Post-processing for the A/e/M unit fix (notes/aem_unit_mismatch.md, step 2).
# "old unit" = the B+ at massfit widths (kpiA_kpi and the fits built on it);
# "matched"  = the B+ at the Z/J/psi fitted-scale widths (kpiC_fitscale{1,5}).
# Idempotent; skips a comparison whose fit is missing rather than failing.
#
# Caution when reading the sigma(A/e/M) plots: they are in PRIOR units, and the
# unit is what changed, so a B+-alone sigma is not comparable between the two
# arms without multiplying by the width. The joint fits share one unit and are.
set -uo pipefail
WREM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$WREM"
CEPH=/ceph/submit/data/user/p/pmlugato/mz
H=$CEPH/histograms; F=$CEPH/fitresults
P=$HOME/public_html/mz/combined_fits_and_mB_prefit_scale_uncs_matched
IMAGE=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest
run() { singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc "cd '$WREM' && source setup.sh >/dev/null 2>&1 && export OPENBLAS_NUM_THREADS=1 && $1" 2>&1 | grep -v "absl\|TensorFlow\|cpu_feature\|rebuild\|FutureWarning\|metadata was not\|alpha stage\|warnings.warn\|^[AeM] unc"; }
have() { [[ -r "$1" ]] || { echo "  skip, missing $(basename $1)"; return 1; }; }
mkdir -p $P/{02_closure,03_bplus_alone,04_mz_joint,05_mB}
# Optional arguments select sections (02 03 04 05); default all.
SECTIONS=${*:-"02 03 04 05"}
want() { [[ " $SECTIONS " == *" $1 "* ]]; }

OLD_B=$F/alllegs_bplus_asimov_fit.hdf5
M1_B=$F/matched1_bplus_asimov_fit.hdf5
NOB=$F/joint_mz_jpsi_bplus_res5_260914_noB_fit.hdf5
OLD_MZ=$F/joint_mz_jpsi_bplus_res5_260914_fit.hdf5
M_MZ=$F/joint_mz_jpsi_bplus_matched_mz_fit.hdf5
JPA=$F/realmc_D_jpsionly_asimov_fit.hdf5; JPD=$F/realmc_D_jpsionly_data_fit.hdf5
OLD_MBA=$F/realmc_C_asimov_fit.hdf5;     OLD_MBD=$F/realmc_A_data_fit.hdf5
M_MBA=$F/joint_mz_jpsi_bplus_matched_mB_noZ_jpsimu_jpsifull_fit.hdf5
M_MBD=$F/joint_mz_jpsi_bplus_matched_mB_noZ_data_jpsimu_jpsifull_fit.hdf5

if want 02; then
echo "=== 02 closure of the B+ histograms, and the template-induced shift per unit"
if have $H/btojpsik_cvhnano_kpiC_fitscale5.hdf5 && have $H/btojpsik_cvhnano_kpiC_fitscale1.hdf5; then
  run "python scripts/btojpsik/aem_unit_closure.py --ref $H/btojpsik_cvhnano_kpiA_kpi.hdf5 \
        --test $H/btojpsik_cvhnano_kpiC_fitscale5.hdf5 --test1 $H/btojpsik_cvhnano_kpiC_fitscale1.hdf5 \
        --outdir $P/02_closure" | tee $P/02_closure/closure.txt
fi
run "python scripts/btojpsik/template_induced_shift.py --params A e M --tensors \
      'old unit, repro chain (res5)=btojpsik_cvhnano_tensor_joint24_res5_260914.hdf5' \
      'matched, repro chain (scale 1/1/1)=btojpsik_cvhnano_tensor_joint24_matched_mz.hdf5' \
      'old unit, mB chain (realmc)=btojpsik_cvhnano_tensor_joint24_mB_realmc_noZ_jpsifull.hdf5' \
      'matched 1/1/1 stretched A,e x5, mB chain=btojpsik_cvhnano_tensor_joint24_matched_mB_noZ_jpsifull.hdf5'" \
  | tee $P/02_closure/template_induced_shift.txt
fi

want 03 && {
echo "=== 03 standalone B+ Asimov"
if have $M1_B; then
  run "python scripts/rabbit/aem_sigma_per_group.py --fits 'old unit (massfit widths)=$OLD_B' \
        'matched, scale 1/1/1=$M1_B' \
        -o $P/03_bplus_alone/aem_sigma_bplus_alone.png --table $P/03_bplus_alone/table_bplus_alone.md \
        --title 'Standalone B+ Asimov, sigma in PRIOR units (the unit differs per arm)'"
fi

}
want 04 && {
echo "=== 04 Z + J/psi + B+ Asimov (repro chain)"
if have $M_MZ; then
  run "python scripts/rabbit/statscale_report.py --fits 'Z+Jpsi (no B)=$NOB' 'old unit=$OLD_MZ' 'matched=$M_MZ' \
        -o $P/04_mz_joint/mz_report.md"
  # statscale_report's --plot is laid out for the statistics-scaling study and
  # draws nothing for these labels; this is the plot for this comparison.
  run "python scripts/rabbit/matched_units_summary_plot.py -o $P/04_mz_joint/mz_summary.png \
        --poi massShiftZ100MeV --poiScale 100 --poiLabel mZ \
        --title 'Z + J/psi + B+ Asimov (repro chain, unit 1e-3 / 1e-2 / 1e-4)' \
        --fit ref asimov 'Z + J/psi' $NOB --fit old asimov 'Z + J/psi + B+' $OLD_MZ \
        --fit matched asimov 'Z + J/psi + B+' $M_MZ"
  run "python scripts/rabbit/aem_sigma_per_group.py --fits 'Z+Jpsi (no B)=$NOB' 'old unit=$OLD_MZ' 'matched=$M_MZ' \
        -o $P/04_mz_joint/aem_sigma_per_group.png --table $P/04_mz_joint/table_aem.md \
        --title 'Z + J/psi + B+ Asimov: what B+ adds to the scale parameters'"
  CFG="--config scripts/rabbit/mz_impacts_plotconfig.py --poi massShiftZ100MeV --scaleImpacts 100 --impactTitle 'Impact on m<sub>Z</sub> [MeV]' -o $P/04_mz_joint --otherExtensions png pdf"
  run "python rabbit/bin/rabbit_plot_pulls_and_impacts.py $M_MZ -r $OLD_MZ --name 'matched' --refName 'old unit' \
        --mode group --grouping mzscale_aem --oneSidedImpacts --showNumbers --noPulls --postfix mz_matched_vs_old $CFG | grep -c png"
fi

}
want 05 && {
echo "=== 05 J/psi + B+ with mB floating"
SPECS=""
for s in "old unit Asimov=$OLD_MBA" "matched Asimov=$M_MBA" "old unit DATA=$OLD_MBD" "matched DATA=$M_MBD"; do
  [[ -r "${s#*=}" ]] && SPECS="$SPECS '$s'"
done
run "python scripts/rabbit/statscale_report.py --poi massShiftB5MeV --poiScale 5 --poiLabel mB --fits $SPECS -o $P/05_mB/mB_report.md"
FITS=""
for s in "ref data|J/psi alone|$JPD" "old asimov|J/psi + B+|$OLD_MBA" "matched asimov|J/psi + B+|$M_MBA" \
         "old data|J/psi + B+|$OLD_MBD" "matched data|J/psi + B+|$M_MBD"; do
  IFS='|' read -r ak lab f <<< "$s"; [[ -r $f ]] && FITS="$FITS --fit $ak '$lab' $f"
done
run "python scripts/rabbit/matched_units_summary_plot.py -o $P/05_mB/mB_summary.png \
      --poi massShiftB5MeV --poiScale 5 --poiLabel mB \
      --title 'J/psi + B+, mB floating (full J/psi chain, unit 5e-3 / 5e-2 / 1e-4)' $FITS"
for arm in asimov data; do
  if [[ $arm == asimov ]]; then M=$M_MBA; O=$OLD_MBA; J=$JPA; else M=$M_MBD; O=$OLD_MBD; J=$JPD; fi
  have $M || continue
  D=$P/05_mB/$arm; mkdir -p $D
  run "python scripts/rabbit/poi_correlation_families.py $M -o $D/corr_mB_families_matched.png --table $D/corr_mB_families_matched.md"
  run "python scripts/rabbit/aem_sigma_per_group.py --poi massShiftB5MeV --poiScale 5 --poiLabel mB \
        --fits 'J/psi alone=$J' 'old unit=$O' 'matched=$M' \
        -o $D/aem_sigma_per_group.png --table $D/table_aem.md \
        --title '$arm, J/psi + B+ with mB floating: what B+ adds to the scale parameters'"
  CFG="--config scripts/rabbit/mz_impacts_plotconfig.py --poi massShiftB5MeV --scaleImpacts 5 --impactTitle 'Impact on m<sub>B</sub> [MeV]' -o $D --otherExtensions png pdf"
  run "python rabbit/bin/rabbit_plot_pulls_and_impacts.py $M -r $O --name 'matched' --refName 'old unit' \
        --mode group --grouping mb_aem --oneSidedImpacts --showNumbers --noPulls --postfix mB_matched_vs_old $CFG | grep -c png"
  run "python rabbit/bin/rabbit_plot_pulls_and_impacts.py $M --mode ungrouped --num 40 --postfix mB_matched_top40 $CFG | grep -c png"
  for reg in 'barrel_eta12 12' 'endcap_eta0 0'; do set -- $reg
    run "python rabbit/bin/rabbit_plot_cov.py $M -o $D --correlation --showNumbers --title '' \
          --config scripts/rabbit/mz_impacts_plotconfig.py --postfix corr4x4_matched_$1 --params massShiftB5MeV A$2 e$2 M$2 >/dev/null"
  done
done
}
echo "=== done"
