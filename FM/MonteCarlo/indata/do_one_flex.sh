#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Run via:"
    echo "MonteCarlo/do_one.csh Sta 3000 8"
    echo "where sta is station with Sta_data dir"
    echo "and 3000 is the number of models to iterate through per jump"
    echo "and 8 is the core number for parallelization"
    exit
fi

sta=$1
nruns=$2
#sta=TEST

#note that the read/write limits to pretty much only 4 threads will ever be used maximum
nthread=$3
export OMP_NUM_THREADS=$nthread
mdir=`pwd`
#mdir='/data/cnliu/hhh/AmbientNoise/JointInversion_Tutorial/MCMC_TestData'
#note that mdir should contain STA_data and MonteCarlo directories, and will contain STA directory after running this script

#be sure that any edits are included by 'make clean' and 'make do_MC_Para' and 'make do_MC_Para_Post_Process_v3' in the dircode directory
#dircode="/data/cnliu/hhh/AmbientNoise/JointInversion_Tutorial/MCMC_Codes/MCMC_weightmisfitDispOnly__gt4pt9_posCrustToMantle_1MsplineLT2Mspline"
dircode=$mdir"/../../Codes/MCMC_weightmisfitDispOnly__gt4pt9_posCrustToMantle_1MsplineLT2Mspline"

# # Recompile the script
# cd $dircode
# make clean
# make do_MC_Para
# make do_MC_Para_Post_process_v3
# cd $mdir

### Input options
ff1=$mdir"/../"${sta}"_data"/${sta}.ph

ff2=$mdir"/../"${sta}"_data"/${sta}.HV

ff4=$mdir"/../"${sta}"_data"/${sta}.RF
gw=2.5

# ff4=$mdir"/in.rf" #if no RF being used, still need to read something in..
# gw=3.0

ff3=$mdir"/../"${sta}"_data"/mod.${sta}

#ff5=$mdir"/"${sta}"_data"/${sta}.gv

echo $ff1
echo $ff2



fcontrol=$mdir"/"${sta}.control


echo model $ff3 > $fcontrol
echo para $mdir"/in.para_"${sta} >> $fcontrol #mohochange

#echo disp 1 1 1 $ff1 >> $fcontrol #1=Rayleigh,1=number of inputs [justphase(no H/V)],1 --> phase designated file
echo disp 1 2 1 $ff1 3 $ff2 >> $fcontrol #1=Rayleigh,2=number of inputs [phase&H/V],1 --> phase file, 3 --> h/v file
#echo disp 1 3 1 $ff1 2 $ff5 3 $ff2 >> $fcontrol #1=Rayleigh,3=number of inputs [phase&group&H/V],1 --> phase file, 2--> group vel, 3 --> h/v file

echo rf $ff4 $gw >> $fcontrol #rf_file gaussian

echo Qmodel $mdir"/"Qmodel.dat >> $fcontrol

# echo p 1 >> $fcontrol #only use dispersion info
echo p 0.5 >> $fcontrol #receiver function and dispersion are used

# echo mono 1 >> $fcontrol #for model not with series of layers
echo mono -1 >> $fcontrol #for model with series of layers? -- unclear if this works

echo monol 2 >> $fcontrol
echo gradl >> $fcontrol

#number of jumps
echo tt 12 >> $fcontrol

#number of models to look at after jump
echo jt $nruns >> $fcontrol
#echo jt 3000 >> $fcontrol

echo outdir $sta $sta  -2 >> $fcontrol
echo n_thread $nthread >> $fcontrol
echo end >> $fcontrol

echo $sta $fcontrol $dircode

if [[ -s $fcontrol ]]; then
echo "control file $fcontrol - ok!"
else
echo "control file $fcontrol is not generated! Stop!"
break;
fi


# Check the output file exit or not.
if ! test -f $output_file_1; then
  # If the file exist!
  echo $output_file_1 ' does not exist! run the 1st time'
  sleep 2
  sh do_one_flex.sh ${sta} 3000 18 2>&1 | tee 'output_'$sta'_tmp.txt'

else
    #Check the report file 1. If the last line have a finish_string_1 then copy to running directory

    # If the output file not exist then run the 1st time
while [[ "$string_1" != "$check_point_1"]] && [["$n_run_time_1" < 100 ]] # run check point only 100 times
do
          ((n_run_time_1++))
          echo "Run the do_MC_Para for the $n_run_time_1 time!"

          time $dircode/do_MC_Para $fcontrol # run the first code
        # cp output_${sta}*txt ${sta}
done
fi

# # 
# # rm $output_file
# # 
# first_check="do_MC_Para done!!!"
# second_check="post process finish!!! ready to plot!"
# charfit1="Hihi"
# charfit2="Hihi"

# n_para_run=0
# # If the output file not exist then run the 1st time
# while [[ "$charfit1" != "$first_check"]] && [["$n_para_run" < 100 ]] # run check point only 100 times
#     do
#         ((n_para_run++))
#         echo "Run the do_MC_Para for the $n_para_run time!"

#     if ! test -f $output_file; then
#         # echo "$output_file not exist! run the 1st time"
#         charfit1="Haha"
#         time $dircode/do_MC_Para $fcontrol # run the first code

#     else

#         charfit1=$( tail -n 1 "$output_file" )
#         echo "charfit1: " $charfit1
#         if [ "$charfit1" == "$first_check" ]; then 
#           break;
#         else
#           echo "string is not match!!!"
#           time $dircode/do_MC_Para $fcontrol # run the first code
#         fi

        
#     fi
    
# done

# # # # Run the do_MC_Para_Post_process_v3
# n_para_run3=0
# if [[ "$charfit1" == "$first_check" ]];then
#     while [[ "$charfit2" != "$second_check"]] && [["$n_para_run3" < 100 ]] # run check point only 100 times
#     do
#         ((n_para_run3++))
#         echo "Run the do_MC_Para for the $n_para_run3 time!"
        
#         if ! test -f output_"$sta".txt; then
#             # echo "output_"$sta".txt not exist! run the 1st time"

#             time $dircode/do_MC_Para_Post_process_v3 $fcontrol > "output_"$sta".txt" ## get posterior information

#         else

#                 charfit2=$( tail -n 1 output_"$sta".txt )
#                 echo $charfit2
#                 if [ "$charfit2" == "$second_check" ]; then 
#                     break;
#                 else
#                     echo "string is not match!!!"
#                     time $dircode/do_MC_Para_Post_process_v3 $fcontrol > "output_"$sta".txt" ## get posterior information 
#                 fi
#                 # time $dircode/do_MC_Para_Post_process_v3 $fcontrol > "output_"$sta".txt" ## get posterior information 
#         fi
#     done
# fi
