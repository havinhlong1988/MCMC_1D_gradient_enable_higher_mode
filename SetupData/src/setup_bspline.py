#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import warnings
import numpy as np
import pandas as pd
import os, sys
import subprocess
import time
import matplotlib.pyplot as plt
#
import parameters
from src import setup_functions
#
# Suppress all warnings
warnings.filterwarnings("ignore")

def setup(stadata,Vel_dir,src_dir,query_sub,data_sub,in_sub):
    # print(stadata)
    # --------------------------------------------------------------------------------------------------------
    # --------------------------------- Data import and preprocessing ----------------------------------------
    # --------------------------------------------------------------------------------------------------------
    # reading the model file:
    if (parameters.vel_mod_uniform_flag==1):
        modfile = os.path.join(Vel_dir,parameters.modfile); 
        # print("<< ",modfile)
    elif (parameters.vel_mod_uniform_flag==0):
        modfile = os.path.join(Vel_dir,parameters.modfile);
        modfile = os.path.join(modfile,f"{stadata['name']}.mod");
        # print("<< ",modfile)
    # Now read the model data 
    if os.path.exists(modfile):
        mod_data = pd.read_csv(modfile,delim_whitespace=True, na_values=0,names=["dep","vs",'vp'])
        mod_data = mod_data.fillna(0); 
    else:
        print("file %s does not exist!",mod_data)
        os.exit(0)
    
    # save out the Vs model file to check
    modfile_out = os.path.join(query_sub,f"{stadata['name']}_depth_vs.txt"); 
    modvs =  mod_data.iloc[:, :2]; 
    #
    modvs['dep'] = modvs['dep'].apply(lambda x: f"{x:5.3f}")
    modvs['vs'] = modvs['vs'].apply(lambda x: f"{x:5.5f}")
    #
    modvs.to_csv(modfile_out,index=False,header=True,sep=" ")
    #
    dVs=[]; dVsdepth=[]
    modstavs = modvs['vs'].to_numpy().astype(np.float32); 
    allgridmoddepth = modvs['dep'].to_numpy().astype(np.float32);
    tmp=allgridmoddepth; 
    # --------------------------------------------------------------------------------------------------------
    # ---------------------------------  in.para_{sta} file produce   ----------------------------------------
    # --------------------------------------------------------------------------------------------------------
    if (parameters.sed_flag==1):
        #5th column
        SedLayerId='0'; 
        CrustLayerId='1'; 
        locals()['inparasedthick'] = np.append(parameters.sediment_thick,[SedLayerId]); 
         
        if (parameters.man_flag==1):
            ManLayerId='2'; 
            for x in range(0,parameters.mantlenpara):
                locals()['inparamantleVs%s' % x] = np.append(parameters.manlte_vel, [ManLayerId, str(x)]); 
                # print('inparamantleVs',x, locals()['inparamantleVs%s' % x])
        #
        
        sed_value = parameters.sed_value; 
        for x in range(0,parameters.sednpara):
            locals()['inparasedVs%s' % x] = np.append(parameters.sediment_vel,[SedLayerId, str(x)])
            # print('inparasedVs',x, locals()['inparasedVs%s' % x])
    else:
        #5th column
        CrustLayerId='0'; 
        if (parameters.man_flag==1):
            ManLayerId='1'; 
            for x in range(0,parameters.mantlenpara):
                locals()['inparamantleVs%s' % x] = np.append(parameters.manlte_vel, [ManLayerId, str(x)]); 
                # print('inparamantleVs',x, locals()['inparamantleVs%s' % x])
        sed_value=modvs['vs'][0];     
        
    for x in range(0,parameters.crustnpara):
        locals()['inparacrustVs%s' % x] = np.append(parameters.crust_vel, [CrustLayerId, str(x)]); 
        # print('inparacrustVs',x, locals()['inparacrustVs%s' % x]) 
        
# ---------------------- Now we find the sediment depth and moho depth ---------------------------
    modstavs[0] = sed_value; 
    for kk in np.arange(1,len(allgridmoddepth)-1):
            dVs.append(abs((modstavs[kk+1]-modstavs[kk-1])/(allgridmoddepth[kk+1]-allgridmoddepth[kk-1])))
            dVsdepth.append((allgridmoddepth[kk+1]+allgridmoddepth[kk-1])*0.5)
    # print(dVsdepth,dVs)

    prominence,prdep=setup_functions.peakdet2(dVs,dVsdepth)
    # print(modstavs)
    dVsmax1index=np.where(prominence[0:4]==min(prominence[0:4]))[0][0]
    # dVsmax2index=(np.where(prominence[4:]==min(prominence[4:]))[0][0])#+5
    # mohoindex=(np.where(np.abs(np.diff(modstavs[1:-1]))>0.5))[0][0]+1
    mohoindex=8
    mohodepth=allgridmoddepth[mohoindex]
    dVsmax1tmp=np.where(dVs[0:4]==max(dVs[0:4]))[0][0]
    # sediment depth and sediment index here
    seddepth,sedindex=setup_functions.find_nearest(allgridmoddepth,dVsdepth[dVsmax1tmp])    
    # mohodepth,mohoindex=setup_functions.find_nearest(allgridmoddepth,prdep[dVsmax2index-1])
    moho_id = mohoindex+1; 
    # Correction the sediment and moho depth based on model settup
    if (parameters.sed_flag==1):
        crustthick=mohodepth-seddepth; 
        sedthick=seddepth; 
    elif (parameters.sed_flag==0):
        sedthick=0; 
        seddepth=0; 
        crustthick=mohodepth; 
    mantlethick=max(allgridmoddepth)-mohodepth
    pertRangCrustThick=np.round((np.float32(parameters.PertRangeCrustThick)*crustthick/100),1)
    # print("pertRangCrustThick: ", pertRangCrustThick)
    if np.isnan(pertRangCrustThick):
        print('error!!!! percpert')
        sys.exit(0); 
    print(">>>seddepth, mohodepth,crustthick, sedindex,mohoindex:",seddepth, mohodepth,crustthick,sedindex,mohoindex)
    parameters.crust_thick[2]=pertRangCrustThick; 
    if (parameters.man_flag==1): # If no manlte then no crustal thickness change.
        locals()['inparacrustthick']=np.append(parameters.crust_thick,CrustLayerId); 
    # ---------------------- in.para_{sta} print out ----------------------------------------------------------------
    # ------------------------------- >>>>>>>> write in.para_{STA} here !!! <<<<<<< ----------------------------------
    fnpara=os.path.join(data_sub, f"in.para_{stadata['name']}")
    print("Write the in.para: ",fnpara)
    # # Remove the file if exist
    # if os.path.exists(fnpara):
    #     os.remove(fnpara); 
    with open(fnpara,'w') as ff0:
        if (parameters.sed_flag==1):
            for x in range(0,parameters.sednpara):
                ff0.write(' '.join(locals()['inparasedVs%s' % x])+'\n')
        for x in range(0,parameters.crustnpara):
            ff0.write(' '.join(locals()['inparacrustVs%s' % x])+'\n')
        if (parameters.man_flag==1):
            for x in range(0,parameters.mantlenpara):
                ff0.write(' '.join(locals()['inparamantleVs%s' % x])+'\n')
        # crustal pertubation part
        if (parameters.sed_flag==1):
            ff0.write(' '.join(locals()['inparasedthick'])+'\n')
        if (parameters.man_flag==1):
            ff0.write(' '.join(locals()['inparacrustthick'])+'\n')
    ff0.close()
    # --------------------------------------------------------------------------------------------------------
    # ---------------------------------  in.data_{sta} file produce   ----------------------------------------
    # --------------------------------------------------------------------------------------------------------   
    fndata = os.path.join(data_sub, f"in.data_{stadata['name']}")
    if (parameters.is_in_data==1):
        # input data name format
        gvfile = os.path.join(data_sub, f"{stadata['name']}.gv")
        phfile = os.path.join(data_sub, f"{stadata['name']}.ph")
        hvfile = os.path.join(data_sub, f"{stadata['name']}.HV")
        rffile = os.path.join(data_sub, f"{stadata['name']}.RF")
        #
        nindata=0
        print("Write the in.data with misfit based on true input data: ",fndata)
        with open(fndata,'w') as ff0:
            if os.path.exists(phfile):
                ff0.write("1 \n")
                nindata+=1; 
                pindata=1; 
            else:
                ff0.write("0 \n")
                pindata=0; 
            #
            if os.path.exists(gvfile):
                ff0.write("1 \n")
                nindata+=1; 
                gindata=1;  
            else:
                ff0.write("0 \n")
                gindata=0; 
            #
            if os.path.exists(hvfile):
                ff0.write("1 \n")
                nindata+=1; 
                hindata=1; 
            else:
                ff0.write("0 \n")
                hindata=0; 
            #
            if os.path.exists(rffile):
                ff0.write("1 \n")
                nindata+=1; 
                rindata=1; 
            else:
                ff0.write("0 \n")
                rindata=0; 
            print(">>>",pindata,gindata,hindata,rindata,nindata)
            if (parameters.is_equal_weight==0):
                print("Write the weighting 0")
                ff0.write("{:d} \n".format(parameters.is_equal_weight))
                #
                ff0.write("{} \n".format(parameters.phw))
                ff0.write("{} \n".format(parameters.gvw))
                ff0.write("{} \n".format(parameters.hvw))
                ff0.write("{} \n".format(parameters.rfw))
            else:
                print("Write the weighting 1")
                ff0.write("{:d} \n".format(parameters.is_equal_weight))
                #
                ff0.write("{:.6f} \n".format(pindata/nindata))
                ff0.write("{:.6f} \n".format(gindata/nindata))
                ff0.write("{:.6f} \n".format(hindata/nindata))
                ff0.write("{:.6f} \n".format(rindata/nindata))

        ff0.close()
    elif (parameters.is_in_data==0):
        print("Write the in.data with optional data misfit style: ",fndata)
        with open(fndata,'w') as ff0:
            ff0.write("{} \n".format(parameters.iph))
            ff0.write("{} \n".format(parameters.igv))
            ff0.write("{} \n".format(parameters.ihv))
            ff0.write("{} \n".format(parameters.irf))
            if (parameters.is_equal_weight==0):
                print("Write the weighting 0")
                ff0.write("{:d} \n".format(parameters.is_equal_weight))
                #
                ff0.write("{} \n".format(parameters.phw))
                ff0.write("{} \n".format(parameters.gvw))
                ff0.write("{} \n".format(parameters.hvw))
                ff0.write("{} \n".format(parameters.rfw))
            else:
                print("Write the weighting 1")
                ff0.write("{:d} \n".format(parameters.is_equal_weight))
                #
                ff0.write("{:.6f} \n".format(pindata/nindata))
                ff0.write("{:.6f} \n".format(gindata/nindata))
                ff0.write("{:.6f} \n".format(hindata/nindata))
                ff0.write("{:.6f} \n".format(rindata/nindata))
        ff0.close()
        if (parameters.igv==1):
            gvfile = os.path.join(data_sub, f"{stadata['name']}.gv")
            if not os.path.exists(gvfile):
                gvfile = os.path.join(in_sub, "in.gv")
                print("---!!!  Warning! Using fake GV data: {}".format(gvfile))
        else:
            gvfile=""; 
        if (parameters.iph==1):
            phfile = os.path.join(data_sub, f"{stadata['name']}.ph")
            if not os.path.exists(phfile):
                phfile = os.path.join(in_sub, "in.ph")
                print("---!!!  Warning! Using fake PH data: {}".format(phfile))
        else:
            phfile=""; 
        if (parameters.ihv==1):
            hvfile = os.path.join(data_sub, f"{stadata['name']}.HV")
            if not os.path.exists(hvfile):
                hvfile = os.path.join(in_sub, "in.hv")
                print("---!!!  Warning! Using fake HV data: {}".format(hvfile))
        else:
            hvfile=""; 
        if (parameters.irf==1):
            rffile = os.path.join(data_sub, f"{stadata['name']}.RF")
            if not os.path.exists(rffile):
                rffile = os.path.join(in_sub, "in.rf")
                print("---!!!  Warning! Using fake RF data: {}".format(rffile))
        else:
            rffile=""; 
    else:
        print("Wrong values of [is_in_data]. Stop!")
        sys.exit(0)
    # --------------------------------------------------------------------------------------------------------
    # -----------------------------------  {sta}.mod file produce   ------------------------------------------
    # --------------------------------------------------------------------------------------------------------

    # --------------------------- Bspline calculation for model ----------------------------------------------
    if (parameters.sed_flag==1):
        vsBspinversion=modstavs[1:moho_id]; depthBspinversion=allgridmoddepth[1:moho_id]
#       mvsBspinversion=modstavs[mohoindex:]; mdepthBspinversion=allgridmoddepth[mohoindex:]
        mvsBspinversion=modstavs[moho_id:]; mdepthBspinversion=allgridmoddepth[moho_id:]
    else:
        vsBspinversion=modstavs[0:moho_id]; depthBspinversion=allgridmoddepth[0:moho_id]
        mvsBspinversion=modstavs[moho_id:]; mdepthBspinversion=allgridmoddepth[moho_id:]
    # run Bspline for crustal
    cmd=src_dir+'/lfBsp_nlay 0 '+str(len(vsBspinversion))+' 2.0 '+str(parameters.crustnpara)+' 4 '+str(len(vsBspinversion))
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    bspfile='B_spline_'
    FCBs=[]
    # maxtrix forming for crustal
    for kk in np.arange(parameters.crustnpara):
        FCBs.append(np.loadtxt(os.path.dirname(src_dir)+'/'+bspfile+str(kk)+'.txt',usecols=(1),unpack=True))
        cmd='mv '+bspfile+str(kk)+'.txt '+bspfile+str(kk)+'_crust.txt'
        pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        pr.wait()
    # run Bspline for mantle
    cmd=src_dir+'/lfBsp_nlay 0 '+str(len(mvsBspinversion))+' 0.75 '+str(parameters.mantlenpara)+' 4 '+str(len(mvsBspinversion))
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    bspfile='B_spline_'
    mFCBs=[]
    # maxtrix forming for mantal   
    for kk in np.arange(parameters.mantlenpara):
        mFCBs.append(np.loadtxt(os.path.dirname(src_dir)+'/'+bspfile+str(kk)+'.txt',usecols=(1),unpack=True))
        cmd='mv '+bspfile+str(kk)+'.txt '+bspfile+str(kk)+'_mantle.txt'
        pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        pr.wait()
        
    # Inverted the coefficient to get the inverted value
    FCGmatrix=np.vstack([FCBs]).T
    mFCGmatrix=np.vstack([mFCBs]).T
    #print(FCGmatrix,mFCGmatrix)
    # upload and invert bsplines
    FCcrustcoeff=np.linalg.lstsq(FCGmatrix,vsBspinversion)[0]
    mFCcrustcoeff=np.linalg.lstsq(mFCGmatrix,mvsBspinversion)[0]
    # solve the equaltion
    #vsfromcrustcoeff=np.dot(Gmatrix,crustcoeff)
    FCvsfromcrustcoeff=np.dot(FCGmatrix,FCcrustcoeff)
    mFCvsfromcrustcoeff=np.dot(mFCGmatrix,mFCcrustcoeff)
    if (parameters.sed_flag==1):
        sedvalue=[]; seddepth=[]
        sedvalue.append(modstavs[0]); seddepth.append(allgridmoddepth[0])
        sedvalue.append(modstavs[1]); seddepth.append(allgridmoddepth[1])
    
        if sedvalue[0]==0:
            sedvalue[0]=0.5
        if sedvalue[1]==0:
            sedvalue[1]=0.
        
# # ----------------- >>>>> More calculation >>>> print to mod.{STA} ---------------------------------
    # sediment
    if (parameters.sed_flag==1):
        sthick=seddepth[len(seddepth)-1]-seddepth[0]
        locals()['modsedlayer']=SedLayerId+' '+parameters.modlinear+' %4.1f'%(sthick)+' %d'%(parameters.sednpara)+' %5.2f'%(float(parameters.dds))+\
        ' %5.6f %5.6f '%(sedvalue[0],sedvalue[1])+parameters.sedrhoflag+' '+parameters.sedQflag+' '+parameters.sedPflag+' '+\
        parameters.sedmodvpvs+' '+parameters.sedVpVs1+' '+parameters.sedDrho+' '+parameters.sedDrho1+' '+parameters.seddVs+' '+\
        parameters.seddvs1+' '+parameters.sedfdvs+' '+parameters.sedfdvs1
    # crust
    printcrusttmp=''
    for kk in np.arange(parameters.crustnpara):
        printcrusttmp=printcrusttmp+' %5.6f'%(FCcrustcoeff[kk])
    # 
    locals()['modcrustlayer']=CrustLayerId+' '+parameters.modsplines+' %4.1f %d %5.2f'%(crustthick,parameters.crustnpara,float(parameters.ddc))+printcrusttmp+\
    ' '+parameters.crustrhoflag+' '+parameters.crustQflag+' '+parameters.crustPflag+' '+parameters.crmodvpvs+' '+parameters.crVpVs1+' '+\
    parameters.crDrho+' '+parameters.crDrho1+' '+parameters.crdVs+' '+parameters.crdvs1+' '+parameters.crfdvs+' '+parameters.crfdvs1
    # Mantle
    if (parameters.man_flag==1):
        printmantletmp=''
        for kk in np.arange(parameters.mantlenpara):
            printmantletmp=printmantletmp+' %5.6f'%(mFCcrustcoeff[kk])
        locals()['modmantlelayer']=ManLayerId+' '+parameters.modsplines+' %4.1f %d %5.2f'%(mantlethick,parameters.mantlenpara,float(parameters.ddm))+printmantletmp+\
        ' '+parameters.mantlerhoflag+' '+parameters.mantleQflag+' '+parameters.mantlePflag+' '+parameters.mtmodvpvs+' '+parameters.mtVpVs1+' '+\
        parameters.mtDrho+' '+parameters.mtDrho1+' '+parameters.mtdVs+' '+parameters.mtdVs1+' '+parameters.mtfdvs+' '+parameters.mtfdvs1
# ------------------------------- >>>>>>>> write mod.{STA} here !!! <<<<<<< ----------------------------------
    fnmod= os.path.join(data_sub, f"mod.{stadata['name']}")
    print("Write the model file in bspline style: ",fnmod)
    with open (fnmod, 'w') as ff2:
        if (parameters.sed_flag==1):
            ff2.write(locals()['modsedlayer']+'\n'); 
        ff2.write(locals()['modcrustlayer']+'\n'); 
        if (parameters.man_flag==1):
            ff2.write(locals()['modmantlelayer']+'\n')
    ff2.close()
# ------------------------------- >>>>>>>> write to in.connector here !!! <<<<<<< ----------------------------------

    fconnect= os.path.join(data_sub, 'in.connector')
    print("Write the connector file to: ",fconnect)
    with open(fconnect,'w') as ff3:
        ff3.write(str(parameters.MC_inversion_type)+"\n")
        ff3.write(str(parameters.RF_gaussian_width)+"\n")
        ff3.write(str(parameters.MC_number_of_jump)+"\n")
        ff3.write(str(parameters.MC_number_of_iteration)+"\n")
        ff3.write(str(parameters.MC_number_of_cores)+"\n")
        ff3.write(str(parameters.MC_percentage_post_process_select)+"\n")
        ff3.write(str(parameters.depth_step)+"\n")
        ff3.write(str(parameters.modout_plot_type)+"\n")
        # ff3.write(str(parameters.sed_mono_check)+"\n")
        # ff3.write(str(parameters.crust_mono_check)+"\n")
    ff3.close()
    # --------------------------------------------------------------------------------------------------------
    # -----------------------------------  Figure plot for initial setup  ------------------------------------
    # --------------------------------------------------------------------------------------------------------
    print("Plot the input model of station %s in Bspline type"%(stadata['name']))
    # read in Vph, H/V, and RF info from file
    if os.path.exists(gvfile):
        allgvper,mergedgv,mergedgvun=np.loadtxt(gvfile,usecols=[0,1,2],unpack=True)
    if os.path.exists(phfile):
        allphper,mergedph,mergedphun=np.loadtxt(phfile,usecols=[0,1,2],unpack=True)
    if os.path.exists(hvfile):
        hvper,hvall,hvun=np.loadtxt(hvfile,usecols=[0,1,2],unpack=True)
    if os.path.exists(rffile):
        rft,rfamp,rfunc=np.loadtxt(rffile,usecols=[0,1,2],unpack=True)
#     ##%% Plot all
# #    '''
    fig=plt.figure(0,figsize=(10,10))
    
    ax1=plt.subplot2grid((3,2), (0,0), rowspan=3) #Velocity models
    ax4=plt.subplot2grid((3,2), (0,1)) #Phase and gv
    ax2=plt.subplot2grid((3,2), (1,1)) #H/V
    ax3=plt.subplot2grid((3,2), (2,1)) #RF
# 
    ax1.set_facecolor('lightyellow')
    ax4.set_facecolor('lightyellow')
    ax2.set_facecolor('lightyellow')
    ax3.set_facecolor('lightyellow')
    
    # Vs
    if (parameters.sed_flag==1):
        ax1.plot(sedvalue,seddepth,'yo',label='Sed Picks')
    ax1.plot(FCvsfromcrustcoeff,depthBspinversion,'ko',label='Est from '+str(parameters.crustnpara)+' Bsplines Vs',alpha=0.4)
    if (parameters.man_flag==1):
        ax1.plot(mFCvsfromcrustcoeff,mdepthBspinversion,'gs',label='Est from '+str(parameters.mantlenpara)+' Bsplines Vs',alpha=0.4)
    ax1.plot(modstavs,allgridmoddepth,'r-',label='Raw Vs Model')
    
    ax1.set_ylim([0,allgridmoddepth[-1]])
    ax1.legend(loc='lower left',fontsize=12)
    
    ax1.set_ylim(ax1.get_ylim()[::-1])
    ax1.set_xlim([0.0,5.7])
    ax1.set_ylabel('Depth (km)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax1.set_xlabel('Vs (km/s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax1.tick_params(labeltop=True)
    ax1.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
    # PHASE
    if os.path.exists(phfile):
        ax2.errorbar(allphper,mergedph,yerr=mergedphun,fmt='.:',color='r',ecolor='r',elinewidth=1.5,capthick=1.5,label='Vph',zorder=5,alpha=0.5)
    if os.path.exists(gvfile):
        ax2.errorbar(allgvper,mergedgv,yerr=mergedgvun,fmt='.:',color='b',ecolor='b',elinewidth=1.5,capthick=1.5,label='Vgv',zorder=5,alpha=0.5)
    # ax2.errorbar(phper,mergedph[0:len(phper)],yerr=mergedphun[0:len(phper)],fmt='.:',color='k',ecolor='k',elinewidth=1.5,capthick=1.5,label='ANT Vph',zorder=5,alpha=0.5)
    #ax2.errorbar(allphper,mergedph,yerr=mergedphun,fmt='.:',color='k',ecolor='r',elinewidth=1.5,capthick=1.5,label='All Vph',zorder=5,alpha=0.5)
    ax2.set_xlabel('Period (s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax2.set_ylabel('Vph (km/s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    #ax2.set_xlim([0.0,105.0])
    ax2.set_ylim([0,5.0])
    ax2.set_xlim([0,30])
    ax2.legend(loc='lower right',fontsize=12)
    ax2.tick_params(labeltop=False)
    ax2.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)

    # H/V
    if os.path.exists(hvfile):
        ax3.errorbar(hvper,hvall,yerr=hvun,fmt='.:',color='r',ecolor='r',elinewidth=1.5,capthick=1.5,label='H/V',zorder=5,alpha=0.5)
    #ax3.errorbar(hvper,hvall,yerr=hvun,fmt='.:',color='r',ecolor='k',elinewidth=1.5,capthick=1.5,label='All H/V',zorder=5,alpha=0.5)
    ax3.set_xlim([0,40])
    ax3.set_xlabel('Period (s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax3.set_ylabel('H/V',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax3.set_ylim([0,2.0])
    ax3.legend(loc='upper right',fontsize=12)
    ax3.tick_params(labeltop=False)
    ax3.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)

    #RF
    if os.path.exists(rffile):
        ax4.errorbar(rft,rfamp,yerr=rfunc,fmt='.:',color='k',ecolor='k',elinewidth=1.5,capthick=0.0,label='RFunc',zorder=5,alpha=0.5)
        ax4.plot(rft,rfamp,'r-',zorder=7,alpha=0.96,label='RF')
        ax4.set_xlim([min(rft)-0.2,max(rft)+0.5])
        ax4.set_ylim([-0.25,0.6])
        ax4.legend(loc='upper right',fontsize=12)
    ax4.set_xlabel('Time (s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax4.set_ylabel('RF Amp',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax4.tick_params(labeltop=False)
    ax4.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)

    plt.suptitle("Starting model and input data for station %s"%(stadata['name']),fontdict = {'family':'serif','color':'darkred','size':25,'weight':'bold'})

    #plt.show()
    figout0 = os.path.join(query_sub,f"StartingModel_{stadata['name']}.png")
    fig.savefig(figout0,bbox_inches='tight',transparent=False,pad_inches=0.1)
    plt.close()
    

    fig, axs = plt.subplots(1,2, figsize=(8,11), sharey=True)
    axs[0].set_facecolor('lightyellow')
    axs[1].set_facecolor('lightyellow')

    axs[1].set_title('B-spline (%s, %s)' %(parameters.crustnpara, parameters.mantlenpara), fontsize=24)

    axs[0].plot(mod_data['vs'], mod_data['dep'], 'ko-',label="Input model")
    if (parameters.sed_flag==1):
        # sed value
        axs[0].plot(sedvalue, tmp[0:2], 'c.-',label="Sed Val")
        # value from initial
        axs[0].plot(vsBspinversion, tmp[1:moho_id], 'ro-',ms=10, alpha=0.5,label="Crustal Val from input")
        if (parameters.man_flag==1):
            axs[0].plot(mvsBspinversion, tmp[moho_id:], 'rs-',ms=10, alpha=0.5,label="Mantle Val from input")
        # value from Bspline coefficient
        axs[0].plot(FCvsfromcrustcoeff, tmp[1:moho_id], 'bo--', alpha=0.5,label="Val from crustal coeff")
        if (parameters.man_flag==1):
            axs[0].plot(mFCvsfromcrustcoeff, tmp[moho_id:], 'bs--', alpha=0.5,label="Val from mantle coeff")
        # value from initial
        axs[1].plot(vsBspinversion, tmp[1:moho_id], 'bo--',ms=10,label="Crustal Val from input")
        if (parameters.man_flag==1):
            axs[1].plot(mvsBspinversion, tmp[moho_id:], 'bs-',ms=10,label="Mantle Val from input")
        # value from Bspline coefficient
        axs[1].plot(FCvsfromcrustcoeff, tmp[1:moho_id], 'ro-',label="Val from crustal coeff")
        if (parameters.man_flag==1):
            axs[1].plot(mFCvsfromcrustcoeff, tmp[moho_id:], 'rs--',label="Val from mantle coeff")

    else: # no sediment
        # value from initial
        axs[0].plot(vsBspinversion, tmp[0:moho_id], 'ro-',ms=10, alpha=0.5,label="Crustal Val from input")
        if (parameters.man_flag==1):
            axs[0].plot(mvsBspinversion, tmp[moho_id:], 'rs-',ms=10, alpha=0.5,label="Mantle Val from input")
        # value from Bspline coefficient
        axs[0].plot(FCvsfromcrustcoeff, tmp[0:moho_id], 'bo--', alpha=0.5,label="Val from crustal coeff")
        if (parameters.man_flag==1):
            axs[0].plot(mFCvsfromcrustcoeff, tmp[moho_id:], 'bs--', alpha=0.5,label="Val from mantle coeff")
        # value from initial
        axs[1].plot(vsBspinversion, tmp[0:moho_id], 'bo-',ms=10,label="Crustal Val from input")
        if (parameters.man_flag==1):
            axs[1].plot(mvsBspinversion, tmp[moho_id:], 'bs--',ms=10,label="Mantle Val from input")
        # value from Bspline coefficient
        axs[1].plot(FCvsfromcrustcoeff, tmp[0:moho_id], 'ro-',label="Val from crustal coeff")
        if (parameters.man_flag==1):
            axs[1].plot(mFCvsfromcrustcoeff, tmp[moho_id:], 'rs--',label="Val from mantle coeff")

    axs[0].tick_params(labeltop=False)
    axs[0].grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
    axs[0].legend(loc='lower left',fontsize=12)
    axs[0].set_ylim(axs[0].get_ylim()[::-1])
    axs[0].set_ylabel('Depth (km)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    axs[0].set_xlabel('Vs (km/s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    
    axs[1].tick_params(labeltop=False)
    axs[1].grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
    axs[1].legend(loc='lower left',fontsize=12)
    axs[1].set_ylim(axs[0].get_ylim())
    axs[1].set_xlim(axs[0].get_xlim())
    axs[1].set_ylabel('Depth (km)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    axs[1].set_xlabel('Vs (km/s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})

    # Adjust layout to prevent overlap
    # plt.tight_layout()
    
    plt.suptitle("Model interpolation for station: %s"%(stadata['name']),fontdict = {'family':'serif','color':'darkred','size':55,'weight':'bold'})

    figout1 = os.path.join(query_sub,f"B-spline_{str(parameters.crustnpara)}_{str(parameters.mantlenpara)}_comparisons_{stadata['name']}.png")
    fig.savefig(figout1,bbox_inches='tight',transparent=False,pad_inches=0.1)
    plt.close()
    
# if __name__ == "__setup_bspline__":
#     setup()
