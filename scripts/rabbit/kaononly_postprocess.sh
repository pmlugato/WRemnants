#!/usr/bin/env bash
# Post-processing for the kaon-only vs all-legs study (notes/jpsi_bplus_mB_fit.md
# option 0). Every output names its arm. Idempotent; skips a comparison whose fit
# is missing rather than failing.
set -uo pipefail
WREM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$WREM"
F=/ceph/submit/data/user/p/pmlugato/mz/fitresults
P=$HOME/public_html/mz/btojpsik_kaononly_260921
IMAGE=/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest
run() { singularity exec --bind /scratch/,/ceph/,/home/,/work/,/tmp/ "$IMAGE" bash -lc "cd '$WREM' && source setup.sh >/dev/null 2>&1 && $1" 2>&1 | grep -v "absl\|TensorFlow\|cpu_feature\|rebuild\|FutureWarning\|metadata was not"; }
have() { [[ -r "$1" ]] || { echo "  skip, missing $(basename $1)"; return 1; }; }

KO_B=$F/kaononly_bplus_asimov_fit.hdf5;      AL_B=$F/alllegs_bplus_asimov_fit.hdf5
KO_MZ=$F/joint_mz_jpsi_bplus_kaononly_mz_fit.hdf5;               AL_MZ=$F/joint_mz_jpsi_bplus_res5_260914_fit.hdf5
KO_MZ10=$F/joint_mz_jpsi_bplus_kaononly_mz_data10_fit.hdf5;      AL_MZ10=$F/joint_mz_jpsi_bplus_statscale_mc2x_data10_fit.hdf5
KO_MZ50=$F/joint_mz_jpsi_bplus_kaononly_mz_data50_fit.hdf5;      AL_MZ50=$F/joint_mz_jpsi_bplus_statscale_mc2x_data50_fit.hdf5
NOB=$F/joint_mz_jpsi_bplus_res5_260914_noB_fit.hdf5
KO_MBA=$F/kaononly_mB_asimov_fit.hdf5;       AL_MBA=$F/realmc_C_asimov_fit.hdf5
KO_MBD=$F/kaononly_mB_data_fit.hdf5;         AL_MBD=$F/realmc_A_data_fit.hdf5
JPD=$F/realmc_D_jpsionly_data_fit.hdf5

echo "=== 01 standalone B+ Asimov: kaononly vs alllegs"
if have $KO_B && have $AL_B; then
  run "python scripts/rabbit/aem_sigma_per_group.py --poi massShiftZ100MeV --fits 'alllegs=$AL_B' 'kaononly=$KO_B' \
        -o $P/01_bplus_alone/aem_sigma_bplus_alone.png --table $P/01_bplus_alone/table_bplus_alone.md \
        --title 'Standalone B+ Asimov: A/e/M variations on all legs vs kaon only'"
fi

echo "=== 02 Z + J/psi + B+ Asimov: kaononly vs alllegs, and the data scalings"
SPECS=""
for s in "Z+Jpsi (no B)=$NOB" "alllegs data1=$AL_MZ" "kaononly data1=$KO_MZ" \
         "alllegs data10=$AL_MZ10" "kaononly data10=$KO_MZ10" \
         "alllegs data50=$AL_MZ50" "kaononly data50=$KO_MZ50"; do
  [[ -r "${s#*=}" ]] && SPECS="$SPECS '$s'"
done
run "python scripts/rabbit/statscale_report.py --fits $SPECS -o $P/02_mz_joint/mz_report.md --plot $P/02_mz_joint/mz_summary.png"
if have $KO_MZ && have $AL_MZ; then
  run "python scripts/rabbit/aem_sigma_per_group.py --fits 'Z+Jpsi (no B)=$NOB' 'alllegs=$AL_MZ' 'kaononly=$KO_MZ' \
        -o $P/02_mz_joint/aem_sigma_per_group.png --table $P/02_mz_joint/table_aem.md \
        --title 'Z + J/psi + B+ Asimov: what B+ adds to the scale parameters'"
  CFG="--config scripts/rabbit/mz_impacts_plotconfig.py --poi massShiftZ100MeV --scaleImpacts 100 --impactTitle 'Impact on m<sub>Z</sub> [MeV]' -o $P/02_mz_joint --otherExtensions png pdf"
  run "python rabbit/bin/rabbit_plot_pulls_and_impacts.py $KO_MZ -r $AL_MZ --name 'kaon-only' --refName 'all legs' \
        --mode group --grouping mzscale_aem --oneSidedImpacts --showNumbers --noPulls --postfix mz_kaononly_vs_alllegs $CFG | grep -c png"
fi

echo "=== 03 J/psi + B+ with mB floating: kaononly vs alllegs, Asimov and data"
SPECS=""
for s in "alllegs Asimov=$AL_MBA" "kaononly Asimov=$KO_MBA" "alllegs DATA=$AL_MBD" "kaononly DATA=$KO_MBD"; do
  [[ -r "${s#*=}" ]] && SPECS="$SPECS '$s'"
done
run "python scripts/rabbit/statscale_report.py --poi massShiftB5MeV --poiScale 5 --poiLabel mB --fits $SPECS -o $P/03_mB_joint/mB_report.md"
if have $KO_MBD && have $AL_MBD; then
  run "python scripts/rabbit/poi_correlation_families.py $KO_MBD -o $P/03_mB_joint/corr_mB_families_kaononly.png --table $P/03_mB_joint/corr_mB_families_kaononly.md"
  run "python scripts/rabbit/aem_sigma_per_group.py --poi massShiftB5MeV --poiScale 5 --poiLabel mB \
        --fits 'Jpsi alone=$JPD' 'alllegs=$AL_MBD' 'kaononly=$KO_MBD' \
        -o $P/03_mB_joint/aem_sigma_per_group.png --table $P/03_mB_joint/table_aem.md \
        --title 'DATA, J/psi + B+ with mB floating: what B+ adds to the scale parameters'"
  CFG="--config scripts/rabbit/mz_impacts_plotconfig.py --poi massShiftB5MeV --scaleImpacts 5 --impactTitle 'Impact on m<sub>B</sub> [MeV]' -o $P/03_mB_joint --otherExtensions png pdf"
  run "python rabbit/bin/rabbit_plot_pulls_and_impacts.py $KO_MBD -r $AL_MBD --name 'kaon-only' --refName 'all legs' \
        --mode group --grouping mb_aem --oneSidedImpacts --showNumbers --noPulls --postfix mB_kaononly_vs_alllegs $CFG | grep -c png"
  run "python rabbit/bin/rabbit_plot_pulls_and_impacts.py $KO_MBD --mode ungrouped --num 40 --postfix mB_kaononly_top40 $CFG | grep -c png"
  for reg in 'barrel_eta12 12' 'endcap_eta0 0'; do set -- $reg
    run "python rabbit/bin/rabbit_plot_cov.py $KO_MBD -o $P/03_mB_joint --correlation --showNumbers --title '' \
          --config scripts/rabbit/mz_impacts_plotconfig.py --postfix corr4x4_kaononly_$1 --params massShiftB5MeV A$2 e$2 M$2 >/dev/null"
  done
fi
echo "=== done; building the index"
python3 scripts/rabbit/kaononly_index.py
