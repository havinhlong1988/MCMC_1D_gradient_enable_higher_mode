#!/bin/bash

# ---------------------------------------------------------------------------
# OpenMP worker-thread stack size.
#
# The Fortran forward model (DISP2/surfa.f) declares large LOCAL arrays -- e.g.
# eight  dimension yy1(nsize,5)..yz4(nsize,5)  plus many nsize=1000 arrays --
# which live on the stack of whichever thread runs them.
#
# Linux gives secondary threads an 8 MB stack, so this always fit on CentOS.
# macOS gives them only 512 KB, so the deeper forward-model paths overflow and
# the run dies with an intermittent "Segmentation fault: 11" / "Bus error: 10"
# (it depends on the data, which is why some stations run fine and others do
# not). Verified: 256K -> SIGBUS immediately, default -> SIGSEGV, 64M -> stable.
#
# Setting this is harmless on Linux (the stack is virtual and committed lazily),
# so the same script works unchanged on CentOS.
# ---------------------------------------------------------------------------
export OMP_STACKSIZE=64M

# -------------------- Run setup the data ---------------------------------------
# cd SetupData/
# python setup_VN.py
# cd ..
# -------------- Recomplile the MCMC code ---------------------------------------
echo "Run the MCMC inversion. Recompile the program or not??? "
echo ""
echo " If yes -> Press y "
echo " If no -> Press other key "
echo ""
# read -p "Press any key to continue... " -n1 -s
echo "Press any key to continue ... "
while [ true ] ; do
# read -t 3 -n 1
read key
if [ "$key" = "y" ] || [ "$key" = "Y" ] ; then
    echo "-------------------------------------"
    echo "Recompile the program !!!"
    echo "-------------------------------------"
    
    cd ./Codes/MCMC_flex
    make clean
    make do_MC_Para
    make do_MC_Para_Post_process_v3
    cd ../../
    break;
else
    echo "-------------------------------------"
    echo "Run without recompile the program !!!"
    echo "-------------------------------------"
    break;
fi
done

# -------------------------------------------------------------------------------
# # Start MCMC inversion

echo "-----------------------------------------------"
echo "Run the MCMC inversion for loop of station !!!"
echo "-----------------------------------------------"

cd CHT/MonteCarlo
num_try=20; # If run false - then try again until $num_try times
for sta in $(awk '{print $1}' ../sta_now)
# for sta in $(awk '{print $1}' ../runstations.lst)
do
# Check if the directory does not exist
if [ ! -d "$sta" ]; then
    echo " > Directory does not exist. Creating directory..."
    # Create the directory
    mkdir -p "$sta"
else
    echo " > Directory $sta exists."
    # echo "Remove old-run directory: ${sta}"
    # rm -r $sta
fi
# Remove the old-output report file
# rm -r output*$sta*


cur_work_dir=`pwd`
# ========================= HVL: Do the MCMC with flexible input: ===========================================
output_file_1='output_'$sta'_tmp.txt'
output_file_2='output_'$sta'.txt'
# check point strings
check_point_1="do_MC_Para done!!!" # do_MC_Para
check_point_2="post process finish!!! ready to plot!" # do_MC_Para_post_process
# pre condition for strings match
string_1="Hehe"
string_2="Haha"
# Number of run 
n_run_time_1=0 # do_MC_Para
n_run_time_2=0 # do_MC_Para_post_process
# ----- MCMC code dir ------
# Main dir:
main_dir=$(dirname $(dirname "$cur_work_dir"))
# Project dir
project_dir=$(dirname "$cur_work_dir")
# Code dir
code_dir=$main_dir/"Codes/MCMC_flex"
# data dir
data_dir=$project_dir"/data/"$sta"_data"
# control link
sfile=$data_dir"/"$sta.control

# >>>>>>>>>>>>>>>>>>>>> make file control <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rm -f $output_file_1
# It important
ls make_control_file* 2>/dev/null
#
echo " > runing the file: 00_make_control_file.sh"
sleep 1
# make the control file for station with 3000 iteration and 20 cores used
bash 00_make_file_control.sh ${sta} 2>&1 | tee $output_file_1
# ################################# >>>>>>>>>>>>>>>>>>>>> runing MCMC part <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
echo "sh make_file_control.sh ${sta} " > 01_run_MCMC.sh
echo "time $code_dir/do_MC_Para $sfile" >> 01_run_MCMC.sh
# Create the report file:
rm -f "00_run_"$sta"_report.txt"
touch "00_run_"$sta"_report.txt"
#
# Start the while true loop
while true; do
  n_run_time_1=$((n_run_time_1 + 1))

  if ! test -f $output_file_1; then
    # If the file exist!
      echo " > $output_file_1 does not exist!"
      touch ${output_file_1}
      # Create for test the loop run
      echo " >> Olala!" >  ${output_file_1}
      #
      sh 01_run_MCMC.sh > "$output_file_1"
    #   sh 01_run_MCMC.sh 2>&1 | tee "$output_file_1"
      # echo " >> Run the do_MC_Para for the first time!"
      echo " >> Run the do_MC_Para for the first time!"  > "00_run_"$sta"_report.txt"
  else
      # If the output file not exist then run the 1st time
      if [ "$n_run_time_1" -le "$num_try" ]; then # run check point only 20 times
          echo " >> $output_file_1 exist! verifying the check point 1"
          echo " >> $output_file_1 exist! verifying the check point 1" >> "00_run_"$sta"_report.txt"
          string_1=$( tail -n 1 $output_file_1 )
          echo $string_1
          sleep 1

          if [ "$string_1" = "$check_point_1" ]; then
              echo " >>> String match, move forward!"
              sleep 1
              echo " >>> String match, move forward!" >> "00_run_"$sta"_report.txt"
              cp $output_file_1 ${sta}/
              break;
          else
              echo " >>> String not match, try again! $string_1 != $check_point_1 "

              echo " >> Run the do_MC_Para for the $n_run_time_1 th time!"
              #
              echo " >> Run the do_MC_Para for the $n_run_time_1 th time!"  >> "00_run_"$sta"_report.txt"

              sh 01_run_MCMC.sh > "$output_file_1"
            #   sh 01_run_MCMC.sh 2>&1 | tee "$output_file_1"
              # ## This one only for check the loop!
              # if [[ "$n_run_time_1" == 10 ]]; then 
              #   echo $check_point_1 > $output_file_1
              # fi
          fi
      fi    
  fi
  # # >>>> Check the report file and copy
  if test -f $output_file_1; then
    cp $output_file_1 ${sta}/
  else
    echo "No such file or directory after run the MCMC_para? $output_file_1"
    echo "No such file or directory after run the MCMC_para? $output_file_1" >> "00_run_"$sta"_report.txt"
    break
  fi
done
# # >>>>>>>>>>>>>>>>>>>>> runing MCMC post process part <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
echo "-----------------------------------------------"
echo "Run the MCMC post-process for station: $sta !!!"
echo "-----------------------------------------------"
# remove the file if have the previous run
rm -f $output_file_2
# >>>>>>>>>>>>>>>
echo "time $code_dir/do_MC_Para_Post_process_v3 $sfile" > 02_run_MCMC_post_process.sh
# Run the post-process.
# NOTE: this binary returns exit code 1 on SUCCESS (codebase convention), so
# `time` prints "Command exited with non-zero status 1" even on success.
# We therefore DO NOT trust the exit code -- we check the finish string instead.
echo " > $output_file_2 does not exist! Running post-process for $sta ..."
sh 02_run_MCMC_post_process.sh > "$output_file_2"
echo "Run the do_MC_Para_post_process_v3 for the first time!" >> "00_run_"$sta"_report.txt"

# >>>> FINISH CONDITION check <<<<
# Search the WHOLE output for the checkpoint string (do NOT use tail -n 1: the
# binary prints "post process finish!!! ready to plot!" followed by a blank line).
if test -f "$output_file_2" && grep -q "$check_point_2" "$output_file_2"; then
    echo " >>> post-process FINISHED ok (checkpoint matched) for $sta"
    echo "post process finish checkpoint matched!" >> "00_run_"$sta"_report.txt"
    cp $output_file_2 ${sta}/
else
    echo " >>> post-process did NOT reach the finish checkpoint! skip plotting $sta"
    echo "post process finish checkpoint NOT found! $output_file_2" >> "00_run_"$sta"_report.txt"
    cp $output_file_2 ${sta}/ 2>/dev/null
    continue
fi

# echo  `ls $sta | wc -l`
# # CNliu fixmantle Weiting HV + PhV
## sh MonteCarlo/do_one_wtVphHV__4pt9kms_fixedMantle_1MsplineLT2Mspline.sh $sta 300 8 > 'output_'$sta'_tmp.txt'

# # ##CNliu Weighting but freemantle + close RF
# sh do_one_wtVphHV__4pt9kms_posCrustToMantle_1MsplineLT2Mspline.sh ${sta} 3000 18 2>&1 | tee 'output_'$sta'_tmp.txt'
# sh do_one_wtVphHV__4pt9kms_posCrustToMantle_1MsplineLT2Mspline.sh ${sta} 3000 4 2>&1 | tee 'output_'$sta'_tmp.txt' # For test
echo "-----------------------------------------------"
echo "Seem all good right now - Plot the figure for station: $sta !!!"
echo "-----------------------------------------------"

# Check the file of post process exist
if ! test -f "$sta/$output_file_2"; then
    echo "$sta/$output_file_2 not exist! Can not plot!"
    echo "$sta/$output_file_2 not exist! Can not plot!" >> "00_run_"$sta"_report.txt"
    # sh do_one_wtVphHV__4pt9kms_posCrustToMantle_1MsplineLT2Mspline.sh ${sta} 3000 18 > 'output_'$sta'_tmp.txt'
else
    python ../../AnalyzeResult/inversion_plot_vfinal_flex.py ${sta} .
fi
# 
    echo "-------------------------------------"
    echo "MCMC process finished for the station $sta !!!"
    echo "-------------------------------------"
done
