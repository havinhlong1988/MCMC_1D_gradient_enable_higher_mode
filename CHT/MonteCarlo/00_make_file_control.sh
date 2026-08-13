#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Run via:"
    echo "MonteCarlo/00_make_file_control Sta"
    echo "where sta is station with Sta_data dir"
    echo "the number of models to iterate through per jump included in "
    echo "8 is the core number for parallelization included in "
    exit
fi

# -----------------------------------------------------------------------------------------------------------------------------------
export OMP_NUM_THREADS=$nthread
curdir=$(pwd)
mdir=$(dirname $(dirname "$curdir"))
datadir="$(dirname "$curdir")/data"
#
sta=$1
# connector file
f0=$datadir"/"${sta}"_data/in.connector"
echo $f0    
#
# gw=`echo ${f0} | awk 'NR==1 { printf "%2.1f" $1}'`
invtype=$(awk 'NR==1 { printf "%.1f\n", $1 }' "$f0") # inversion style 1  or 2
gw=$(awk 'NR==2 { printf "%.1f\n", $1 }' "$f0") # gaussian width for RF
njumps=$(awk 'NR==3 { printf "%d\n", $1 }' "$f0") # Number of jump
nruns=$(awk 'NR==4 { printf "%d\n", $1 }' "$f0") # number of inversion model
nthread=$(awk 'NR==5 { printf "%d\n", $1}' "$f0") # number of cores use
selstyle=$(awk 'NR==6 { printf "%d\n", $1}' "$f0") # selection style: -1 percentage, 1 absolute(min misfit + value)
mc_range=$(awk 'NR==7 { printf "%.2f\n", $1}' "$f0") # selection threshold value (percent if selstyle=-1, absolute offset if selstyle=1)
dstep=$(awk 'NR==8 { printf "%.1f\n", $1}' "$f0") # average step for final model
plottype=$(awk 'NR==9 { printf "%d\n", $1}' "$f0") # plot style layercake or gradient
mc_qc=$(awk 'NR==10 { printf "%d\n", $1}' "$f0") # model quality control flag (use goodmodel function to constrain the model)
sedmonocheck=$(awk 'NR==11 { printf "%d\n", $1}' "$f0") # check the monochromatic increment of sedimentary
crustmonocheck=$(awk 'NR==12 { printf "%d\n", $1}' "$f0") # check the monochromatic increment of crust

# -----------------------------------------------------------------------------------------------------------------------------------

#note that curdir should contain STA_data and MonteCarlo directories, and will contain STA directory after running this script
codedir=$mdir"/Codes/MCMC_flex"
# now print the line to find the values
### Input options
ff0=$datadir"/"${sta}"_data"/${sta}.ph
ff1=$datadir"/"${sta}"_data"/${sta}.gv
ff2=$datadir"/"${sta}"_data"/${sta}.HV
ff3=$datadir"/"${sta}"_data"/${sta}.RF
ff4=$datadir"/"${sta}"_data/mod."${sta}
ff5=$datadir"/"${sta}"_data/in.data_"${sta} # data avaiable flag file
ff6=$curdir"/indata/Qmodel.dat" # Qmodel
# -----------------------------------------------------------------------------------------------------------------------------------
fcontrol=$datadir"/"$sta"_data/"${sta}.control
touch $fcontrol

echo model $ff4 > $fcontrol
echo para $datadir"/"$sta"_data/in.para_"${sta} >> $fcontrol #mohochange
# Check the in.data_{sta} file
phf=$(head -n 1 "$ff5" | tr -s ' ' | sed -E 's/^ | $//g')
gvf=$(awk 'NR==2' "$ff5" | tr -s ' ' | sed -E 's/^ | $//g')
hvf=$(awk 'NR==3' "$ff5" | tr -s ' ' | sed -E 's/^ | $//g')
rff=$(awk 'NR==4' "$ff5" | tr -s ' ' | sed -E 's/^ | $//g')
# Now based on the data avaiable file to enter the file
echo $phf $gvf $hvf $rff
# set the phase velocity file name
if [ "$phf" -eq 1 ]; then
    echo "Phase file exist!"
    ff0=$ff0
else
    ff0=$curdir"/indata/in.ph"
fi
# set the group velocity file name
if [ "$gvf" -eq 1 ]; then
    echo "Group file exist!"
    ff1=$ff1
else
    ff1=$curdir"/indata/in.gv"
fi
# set the ellipticity file name
if [ "$hvf" -eq 1 ]; then
    echo "HV file exist!"
    ff2=$ff2
else
    ff2=$curdir"/indata/in.hv"
fi
# set the receiver function file name
if [ "$rff" -eq 1 ]; then
    echo "RF file exist!"
    ff3=$ff3
    pval=0.5 # use receiver function and dispersion
else
    ff3=$curdir"/indata/in.rf"
    pval=1 # use dispersion only
fi
# echo disp 1 3 1 $ff0 2 $ff1 3 $ff2 >> $fcontrol #1=Rayleigh,3=number of inputs [phase&group&H/V],1 --> phase file, 2--> group vel, 3 --> h/v file
echo disp 1 2 1 $ff0 3 $ff2 >> $fcontrol #1=Rayleigh,3=number of inputs [phase&group&H/V],1 --> phase file, 2--> group vel, 3 --> h/v file
echo rf $ff3 $gw >> $fcontrol #rf_file gaussian
echo Qmodel $ff6 >> $fcontrol
# echo p 1 >> $fcontrol #only use dispersion info
echo p $pval >> $fcontrol #receiver function and dispersion are used
#
if [ "$mc_qc" -eq 1 ]; then
	echo "Use CALmodel/goodmodel function to constrain model behavior"
    echo mono 1 >> $fcontrol #for model not with series of layers
else
	echo "Model are freely walk - if you want to change, check SetupData/parameter.py - MC_QC_thres"
    echo mono -1 >> $fcontrol #for model with series of layers? -- unclear if this works
fi
# 
echo monol 2 >> $fcontrol
echo gradl >> $fcontrol
#number of jumps
echo tt $njumps >> $fcontrol
#number of models to look at after jump
echo jt $nruns >> $fcontrol
#echo jt 3000 >> $fcontrol
echo outdir $sta $sta  -2 >> $fcontrol
echo n_thread $nthread >> $fcontrol
echo indata $ff5 >> $fcontrol
echo selstyle $selstyle >> $fcontrol
echo mc_range $mc_range >> $fcontrol
echo dstep $dstep >> $fcontrol
echo plottype $plottype >> $fcontrol
echo invtype $invtype >> $fcontrol
echo end >> $fcontrol
# echo "station: " $sta 
# echo "control file location: " $fcontrol 
# echo "code directory: " $codedir
# echo "para file: " $ff5
# ========================================================
if [ -s "$fcontrol" ]; then
echo "control file $fcontrol - ok!"
else
echo "control file $fcontrol is not generated! Stop!"
exit 1
fi