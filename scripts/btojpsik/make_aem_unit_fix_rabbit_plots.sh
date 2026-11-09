#!/bin/bash
# rabbit plots for the A/e/M unit-fix deck; run inside the container (source setup.sh).
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
F=/ceph/submit/data/user/p/pmlugato/mz/fitresults
O=slides/aem_unit_fix/assets/rabbit; mkdir -p $O
OLD_D=$F/realmc_A_data_fit.hdf5;  NEW_D=$F/joint_mz_jpsi_bplus_matched_mB_noZ_data_jpsimu_jpsifull_fit.hdf5
OLD_A=$F/realmc_C_asimov_fit.hdf5; NEW_A=$F/joint_mz_jpsi_bplus_matched_mB_noZ_jpsimu_jpsifull_fit.hdf5
OLD_Z=$F/joint_mz_jpsi_bplus_res5_260914_fit.hdf5; NEW_Z=$F/joint_mz_jpsi_bplus_matched_mz_fit.hdf5
CFG=scripts/rabbit/mz_impacts_plotconfig.py
PI="python rabbit/bin/rabbit_plot_pulls_and_impacts.py"
B="--config $CFG --poi massShiftB5MeV --scaleImpacts 5 --impactTitle 'Impact on m<sub>B</sub> [MeV]' -o $O --otherExtensions png pdf"
eval $PI $NEW_D -r $OLD_D --name "'matched unit'" --refName "'old unit'" --mode group --grouping mb_aem --oneSidedImpacts --showNumbers --noPulls --postfix mB_data $B
eval $PI $NEW_A -r $OLD_A --name "'matched unit'" --refName "'old unit'" --mode group --grouping mb_aem --oneSidedImpacts --showNumbers --noPulls --postfix mB_asimov $B
eval $PI $NEW_Z -r $OLD_Z --name "'matched unit'" --refName "'old unit'" --mode group --grouping mzscale_aem --oneSidedImpacts --showNumbers --noPulls --postfix mZ_asimov \
     --config $CFG --poi massShiftZ100MeV --scaleImpacts 100 --impactTitle "'Impact on m<sub>Z</sub> [MeV]'" -o $O --otherExtensions png pdf
for p in e A; do
  eval $PI $NEW_D -r $OLD_D --name "'matched unit'" --refName "'old unit'" --mode ungrouped --filters "'^${p}[0-9]+\$'" -s label -n 24 --oneSidedImpacts --postfix ${p}_pulls_data $B
done
for f in old:$OLD_D matched:$NEW_D; do
  python rabbit/bin/rabbit_plot_cov.py ${f#*:} -o $O --correlation --showNumbers --title '' --config $CFG \
      --postfix corr4x4_${f%%:*}_eta12 --params massShiftB5MeV A12 e12 M12 >/dev/null
done
ls $O | grep png
# barrel-only e_i pulls (eta groups 6-17) for the main slide: 12 rows display larger than 24
eval $PI $NEW_D -r $OLD_D --name "'matched unit'" --refName "'old unit'" --mode ungrouped --filters "'^e([6-9]|1[0-7])\$'" -s label -n 12 --oneSidedImpacts --postfix e_pulls_barrel_data $B
