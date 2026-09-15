# Project directory: pass it as the first argument, e.g. "bash 05_do_replot_MCMC.sh CHT".
# Defaults to FM (the active project); pass CHT to go back to the old one.
project_dir="${1:-FM}"
cd "$project_dir/MonteCarlo"

for sta in $(awk '{print $1}' ../sta_now)
do
output_file_2='output_'$sta'.txt'
if ! test -f "$sta/$output_file_2"; then
    echo "$sta/$output_file_2 not exist! Can not plot!"
    echo "$sta/$output_file_2 not exist! Can not plot!" >> "00_run_"$sta"_report.txt"
    # sh do_one_wtVphHV__4pt9kms_posCrustToMantle_1MsplineLT2Mspline.sh ${sta} 3000 18 > 'output_'$sta'_tmp.txt'
else
    python ../../AnalyzeResult/inversion_plot_vfinal_flex.py ${sta} .
    # zoom-in Vph fit (fundamental + higher mode) with posterior background
    python ../../AnalyzeResult/plot_hpmode_check.py ${sta} . || true
fi

done
