#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
# --------------------------------- Parameter to settup the run for MCMC inversion can be here -----------------------------------
# print("import success!")
##################################################################################################################################
#
# Structure need in each Project_direcotry:
# + data: The data all plan here. The whole setup data - MCMC run and plot output will be navigate to those directory
# + Link to data: data_dir = {Project_dir}/data/{sta}_data
#       * {data_dir}/{sta}.gv : group velocity dispersion file
#       * {data_dir}/{sta}.ph : phase velocity dispersion file
#       * {data_dir}/{sta}.HV : elipticity dispersion file
#       * {data_dir}/{sta}.RF : receiver function file
#   Note: As least 1 type of data need to exist!
##################################################################################################################################
### ==============================================================================================================================
##                    mod.{sta}
### ==============================================================================================================================
# ------------------ seting for the model type  ------------------------------------------
# mod_type_flag: 1 = layered; 2 = Bspline ;
# Note: if mode_type_flag = 2 and sediment = 1 then sediment model type = linear
mod_type_flag = 1
# ---------------------------------------------------------------------------------------
# model flag when use similar 1D velocity model or difference velocity layer for each survey points
# (vel_mod_uniform_flag = 0: use difference models - need define model names consistent with survey point - see modfile_dir)
# (vel_mod_uniform_flag = 1: use 1 model for all)
vel_mod_uniform_flag = 1
# sed_flag: = 0 no sediment setting, = 1 with sediment setting
# if  sed_flag = 1: the 1st layer of 1D model will be replace by sed_value
sed_flag = 0
# if  man_flag = 0: no mantle setup for run (modify mod STA.mod and inpara.STA files)
# If man_flag = 0 then no depth pertubation for crustal
man_flag = 1
# -----------------------------------------------------------------------------------
#general settings
sed_value=2.0 # fixed sediment value [if sed_flag =0 then no use] 
# number of b-splines in each layer use if mod_type_flag = 2 (leave as an int)
# number of sediment parameter - 4th col in mod.STA (required as assuming a top linear layer here)
if mod_type_flag == 2:
  sednpara=2; # linear
  crustnpara=6; # 6-paras spline
  mantlenpara=4; # 4-paras spline
# Set the variation of the sublayers thickness (total thickness ratio remain ==1)
# yes = 1; no = 0;
sublay_thickchange=1; 
if sublay_thickchange==1:
    '''
    Define the variation of the thickness ratio sublayer for sediment - crust and mantle
    1 for yes and 0 for no
    Current settup for CHT array, only crustall thickness are variation
    '''
    # -------------------------------------
    subsedthickchange=0; 
    subcrustthickchange=1; 
    submantlethickchange=0; 
####################################################################################################################################
# the station name and station coordinantes (file name only):
stafilename='station_cor.lst'
# ------------- velocity model filename setting here (file name only) -------------------------------
# if vel_mod_uniform_flag = 1, use the 1D model define here!
if (vel_mod_uniform_flag==1):
    # modfile='NTW1d_H14_modified.dat'
    modfile='NTW_1D_Liu21_modified.dat'
    # modfile='phase1_VsV.dat'
    # modfile='../Vel_mod/NTW1d_H14_1km'
elif (vel_mod_uniform_flag==0): # set directory name only
    modfile='vel_mods_CHT'
else:
    print("Wrong vel_mod_uniform_flag, stop!")
    sys.exit()
####################################################################################################################################
### ==============================================================================================================================
##                    in.data_{sta} option
### ==============================================================================================================================
# Control the in.data_{sta} file, these are two options.
# + 1: Based on the data avaiable to calculate the misfit (e.g. you have Vph and RF only, then misfit based on them)
# + 0: Based on data you want to fit (e.g. you have full 4 types data Vph, Vgv, HV and RF but want to use 1 or coulple of the for misfit)
is_in_data=1; # use the data avaiable
is_equal_weight=0 # use equal weight

if (is_in_data==0):
    '''
    Beware of this option:
    + If you dont have the data each of them, set to 0. If set to 1, it will use fake data at MonteCarlo/indata
    + If you dont want to use the data in misfit calculation -> set to 0;
    '''
    iph=1; # Phase velocity
    igv=0; # Group velocity
    ihv=1; # HV data
    irf=0; # Receiver function data
if (is_equal_weight==0):
    '''
    set different weight for each data type.
    + If data not avaiable then weigh = 0
    + Total weight = 1 (e.g 0.5 0.25 0 0.25)
    '''
    phw=1.0
    gvw=0.0
    hvw=0.0
    rfw=0.0

####################################################################################################################################
### ==============================================================================================================================
##                    in.para_{sta} (Should change here!)
### ==============================================================================================================================
# 1st column options (fix thickness and pertube velocity = 0 | fix velocity and pertube thick  = 1)
PertTypeSedVel='0'; PertTypeSedThick='1'; PertTypeSedSubThick='1'; 
PertTypeCrustVel='0'; PertTypeCrustThick='1'; PertTypeCrustSubThick='1'; 
PertTypeMantleVel='0'; PertTypeMantleThick='1'; PertTypeMantleSubThick='1'; 
# 2nd col options (percent = -1 vs absolute = 1)
PertStyleSedVel='1'; PertStyleSedThick='-1';  PertStyleSedSubThick='-1';  
PertStyleCrustVel='-1'; PertStyleCrustThick='1'; PertStyleCrustSubThick='-1'; 
PertStyleMantleVel='-1';  PertStypeMantleThick='-1'; PertStyleMantleSubThick='-1'; 
# 3rd col; perturbation range
PertRangeSedVel='2.0'; PertRangeSedThick='100'; PertRangeSedSubThick='100'; 
PertRangeCrustVel='90'; PertRangeCrustThick='60'; PertRangeCrustSubThick='90'; 
PertRangeMantleVel='50'; PertRangeMantleSubThick='90'; 
# 4th col #gaussian step width
gwStepSedVel='0.05'; gwStepSedThick='0.1'; gwStepSedSubThick='0.05'; 
gwStepCrustVel='0.05'; gwStepCrustThick='1.0'; gwStepCrustSubThick='0.05'; 
gwStepMantleVel='0.05'; gwStepMantleSubThick='0.05';   
# Notice: The sublayers variation based on the ratio to over sed-crust-mantle total thickness. So variation step (gwStep) is smaller
### =================================== No need to change anything here ========================== 
#### Now setup for in.para_{sta} based above values (colum 5 will set by other script)
sediment_vel = [PertTypeSedVel,PertStyleSedVel,PertRangeSedVel,gwStepSedVel]
#
sediment_thick = [PertTypeSedThick,PertStyleSedThick,PertRangeSedThick,gwStepSedThick]
# 
sediment_sub_thick = [PertTypeSedSubThick,PertStyleSedSubThick,PertRangeSedSubThick,gwStepSedSubThick]
#
crust_vel = [PertTypeCrustVel,PertStyleCrustVel,PertRangeCrustVel,gwStepCrustVel]
# leave the value as percentage here. Will change to absolutely later
crust_thick = [PertTypeCrustThick,PertStyleCrustThick,PertRangeCrustThick,gwStepCrustThick]
#
crust_sub_thick = [PertTypeCrustSubThick,PertStyleCrustSubThick,PertRangeCrustSubThick,gwStepCrustSubThick]
#
manlte_vel = [PertTypeMantleVel,PertStyleMantleVel,PertRangeMantleVel,gwStepMantleVel]
#
mantle_sub_thick = [PertTypeMantleSubThick,PertStyleMantleSubThick,PertRangeMantleSubThick,gwStepMantleSubThick]

### ==============================================================================================================================
##                    mod.{sta} (Should change here!) - See E.M.B's tutorial for more information
### ==============================================================================================================================
modlay="1"; modsplines='2'; modlinear='4'; #2nd col Keep this - do not change
#
sedrhoflag='1'; sedQflag='2'; sedPflag='1'; sedmodvpvs='0'; sedVpVs1='0'; #cols -11, -10, -9 for sediment
crustrhoflag='1'; crustQflag='3'; crustPflag='3'; crmodvpvs='1.70'; crVpVs1='1.73'; #cols -11, -10, -9 for crustal
mantlerhoflag='2'; mantleQflag='4'; mantlePflag='3'; mtmodvpvs='1.75'; mtVpVs1='1.78'; #cols -11, -10, -9 for mantle
#leave these as 0's for now; untested as to if they carry through codes correctly or not
sedDrho='0'; sedDrho1='0'; seddVs='0'; seddvs1='0'; sedfdvs='0'; sedfdvs1='0'
crDrho='0'; crDrho1='0'; crdVs='0'; crdvs1='0'; crfdvs='0'; crfdvs1='0'
mtDrho='0'; mtDrho1='0'; mtdVs='0'; mtdVs1='0'; mtfdvs='0'; mtfdvs1='0'

# Depth step for input model (sediment, crust, mantle)
# Affect only layer-cake and linear model group (Beware, now the forward calculation can dealing with 100 sublayer only)
# (I've test with higher code for 100, 200, 250, 500 and 1000 sublayer, but 100 give the best)
dds=0.2; # sub-layer thickness for sediment (e.g. 0.1 km each)
ddc=0.05; # sub-layer thickness for sediment (e.g. 0.5 km each)
ddm=2.5; # sub-layer thickness for sediment (e.g. 5.0 km each)
### ==============================================================================================================================
##            {sta}.control file parameters: see file MonteCarlo/00_make_file_control.sh for more information
##            Generate when run the 02_do_MCMC.sh 
##            Temporary print to in.connector, then print to {sta}.control by 00_make_file_control.sh
### ==============================================================================================================================
MC_inversion_type = 1; # 1 for layered, 2 for bspline
RF_gaussian_width=2.5; # Gaussian width for foward RFs (need to consistent with factor use to calculate input RFs)
MC_number_of_jump=12; # Number of jump
MC_number_of_iteration=2000; # Number of iteration each run
MC_number_of_cores=5; # Number of cores use in distribution calculation
# Posterior selection style:
#   -1 = by percentage of posterior population (MC_post_process_select in %, e.g. 10 = best 10%)
#    1 = by absolute value: keep models with misfit <= (minimum misfit + MC_post_process_select)
selection_style = -1;
# Model selection threshold value.
#   if selection_style == -1 : percent for all passed MCMC_para (e.g. 20 = best 20%, max 100)
#   if selection_style ==  1 : absolute misfit offset above the minimum misfit (e.g. 0.5)
MC_post_process_select=10;
# Final model depth step characterize by average whole selected model.
depth_step=0.05; 
# Output model plot stype (layercake style = 0; gradient style = 1)
modout_plot_type=1; 
# Threshold to use the goodmodel constrain function (1) or not (0) (force model in some way. Kinly check CALmodel.C/goodmodel) 
MC_QC_thres=1; 
############################################################ END ########################################################################
