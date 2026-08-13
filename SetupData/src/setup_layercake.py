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
from src import setup_functions, setup_bspline, setup_layercake, setup_gradient
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
        print(f"file {modfile} does not exist!")
        sys.exit(0)
    
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
    # ---------------------- Now we find the sediment depth and moho depth ---------------------------
    if (parameters.sed_flag==1):
        modstavs[0] = parameters.sed_value;  
    for kk in np.arange(1,len(allgridmoddepth)-1):
            dVs.append(abs((modstavs[kk+1]-modstavs[kk-1])/(allgridmoddepth[kk+1]-allgridmoddepth[kk-1])))
            dVsdepth.append((allgridmoddepth[kk+1]+allgridmoddepth[kk-1])*0.5)
    # print(dVsdepth,dVs)

    prominence,prdep=setup_functions.peakdet2(dVs,dVsdepth)
    print(prominence,prdep)
    dVsmax1index = (np.where(prominence[:2] == np.array(prominence[:2]).min())[0][0] if len(prominence[:2]) else 0)
    dVsmax1index = (np.where(np.asarray(prominence)[:2] == min(prominence[:2]))[0][0] if len(prominence[:2]) else 0)
    dVsmax1tmp   = (np.where(np.asarray(dVs)[:2] == np.asarray(dVs)[:2].max())[0][0] if len(dVs[:2]) else 0)

    # sediment depth and sediment index here
    seddepth,sedindex=setup_functions.find_nearest(allgridmoddepth,dVsdepth[dVsmax1tmp])    
    # mohodepth,mohoindex=setup_functions.find_nearest(allgridmoddepth,prdep[dVsmax2index])
    print(allgridmoddepth,modstavs)
    # time.sleep(1)
    # mohodepth, mohoindex = setup_functions.detect_moho_from_layer_cake(
    if (parameters.man_flag==1):
        mohodepth, mohoindex = setup_functions.detect_moho_duplicate_depth(
        allgridmoddepth,
        modstavs,
        min_jump=0.2
        )
        print(mohodepth,mohoindex)
        # time.sleep(1)
        moho_id = mohoindex+1; 
    else:
        mohodepth = allgridmoddepth.max()
        mohoindex = int(np.argmax(allgridmoddepth))
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
    print(">>>seddepth, mohodepth,crustthick, sedindex,mohoindex:",seddepth, mohodepth,crustthick,sedindex,mohoindex)
    # ---------------------------------------------------------------------------------------------------------------
    # we calculate the sednpara, crustnpara, mantlenpara based on sedindex and mohoindex above:     
    if (parameters.sed_flag==1): # If have sediment setting
        sednpara = sedindex; 
        if sednpara==0: sednpara=1; 
        crustnpara = (mohoindex - sedindex); 
        mantlenpara = (len(allgridmoddepth)-mohoindex)-2; 
        if (parameters.man_flag==0):
            mantlenpara = 0
    elif (parameters.sed_flag==0): # If no sediment setting
        sednpara=0
        crustnpara = mohoindex#+1;
        mantlenpara = (len(allgridmoddepth)-mohoindex)-2; 
        if (parameters.man_flag==0):
            mantlenpara = 0
    else:
        print("Wrong set_flag setting. Stop!",parameters.set_flag)
        sys.exit(0)
    print("check n-para sed, crust, mantle: ", sednpara, crustnpara, mantlenpara)
    # ---------------------------------------------------------------------------------------------------------------
    
    if (parameters.sed_flag==1):
        #5th column
        SedLayerId='0'; 
        CrustLayerId='1'; 
        locals()['inparasedthick'] = np.append(parameters.sediment_thick,[SedLayerId]); 
        
        if (parameters.man_flag==1):
            ManLayerId='2'; 
            for x in range(0,mantlenpara):
                locals()['inparamantleVs%s' % x] = np.append(parameters.manlte_vel, [ManLayerId, str(x)]); 
            
                if parameters.sublay_thickchange == 1: # If allow the sublayer thickness variation 
                    if parameters.submantlethickchange == 1:
                        locals()['inparamantlesubthick%s' % x] = np.append(parameters.mantle_sub_thick,[ManLayerId, str(x)]); 
            #
        
        for x in range(0,sednpara):
            locals()['inparasedVs%s' % x] = np.append(parameters.sediment_vel,[SedLayerId, str(x)])
            # print('inparasedVs',x, locals()['inparasedVs%s' % x])
            if parameters.sublay_thickchange == 1: # If allow the sublayer thickness variation 
                if parameters.subsedthickchange == 1:
                    locals()['inparasedsubthick%s' % x] = np.append(parameters.sediment_sub_thick,[SedLayerId, str(x)]); 
    else:
        #5th column
        CrustLayerId='0'; 
        if (parameters.man_flag==1):
            ManLayerId='1'; 
            for x in range(0,mantlenpara):
                locals()['inparamantleVs%s' % x] = np.append(parameters.manlte_vel, [ManLayerId, str(x)]); 
                # print('inparamantleVs',x, locals()['inparamantleVs%s' % x])   
                if parameters.sublay_thickchange == 1: # If allow the sublayer thickness variation 
                    if parameters.submantlethickchange == 1:
                        locals()['inparamantlesubthick%s' % x] = np.append(parameters.mantle_sub_thick,[ManLayerId, str(x)]);  
        
    for x in range(0,crustnpara):
        locals()['inparacrustVs%s' % x] = np.append(parameters.crust_vel, [CrustLayerId, str(x)]); 
        # print('inparacrustVs',x, locals()['inparacrustVs%s' % x]) 
        if parameters.sublay_thickchange == 1: # If allow the sublayer thickness variation 
            if parameters.subcrustthickchange == 1:
                locals()['inparacrustsubthick%s' % x] = np.append(parameters.crust_sub_thick,[CrustLayerId, str(x)]); 
        
    # Crust thickness calculation
    pertRangCrustThick=np.round((np.float(parameters.PertRangeCrustThick)*crustthick/100),2)
    # print("pertRangCrustThick: ", pertRangCrustThick)
    if np.isnan(pertRangCrustThick):
        print('error!!!! percpert to pertRangCrustThick')
        sys.exit(0); 

    parameters.crust_thick[2]=pertRangCrustThick; 
    if (parameters.man_flag==1): # If no manlte then no crustal thickness change.
        locals()['inparacrustthick']=np.append(parameters.crust_thick,CrustLayerId); 
    
        

    
    # ------------------------------- >>>>>>>> write in.para_{STA} here !!! <<<<<<< ----------------------------------
    fnpara=os.path.join(data_sub, f"in.para_{stadata['name']}")
    print("Write the in.para: ",fnpara)
    # # Remove the file if exist
    # if os.path.exists(fnpara):
    #     os.remove(fnpara); 
    with open(fnpara,'w') as ff0:
        if (parameters.sed_flag==1):
            for x in range(0,sednpara):
                ff0.write(' '.join(locals()['inparasedVs%s' % x])+'\n')
        for x in range(0,crustnpara):
            ff0.write(' '.join(locals()['inparacrustVs%s' % x])+'\n')
        if (parameters.man_flag==1):
            for x in range(0,mantlenpara):
                ff0.write(' '.join(locals()['inparamantleVs%s' % x])+'\n')
        # crustal pertubation part
        if (parameters.sed_flag==1):
            ff0.write(' '.join(locals()['inparasedthick'])+'\n')
        if (parameters.man_flag==1):
            ff0.write(' '.join(locals()['inparacrustthick'])+'\n')
        if parameters.sublay_thickchange == 1:
            if (parameters.sed_flag==1):
                if parameters.subsedthickchange == 1:
                    for x in range(0,sednpara):
                        ff0.write(' '.join(locals()['inparasedsubthick%s' % x])+'\n')
            if parameters.subcrustthickchange == 1:
                for x in range(0,crustnpara):
                    ff0.write(' '.join(locals()['inparacrustsubthick%s' % x])+'\n')
            if (parameters.man_flag==1):
                if parameters.submantlethickchange == 1:        
                    for x in range(0,mantlenpara):
                        ff0.write(' '.join(locals()['inparamantlesubthick%s' % x])+'\n')
                    

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

    # ------------------------ calculate through the model layers -----------------------------------------------
    zero_d = allgridmoddepth[0]; 
    # 
    tmpsv0=''
    tmpsd0=''
    # 
    tmpcv0=''
    tmpcd0=''
    # 
    tmpmv0=''
    tmpmd0=''
    # forming layered model or linear model
    mod_vs = []
    mod_d = []
    for i,depth in enumerate(allgridmoddepth):   
        # forming for plot and mod.STA
        if (i==0):
            if (parameters.sed_flag==1):
                mod_d.append(depth)
                mod_vs.append(parameters.sed_value)
            elif (parameters.sed_flag==0):
                mod_d.append(depth)
                mod_vs.append(float(modvs['vs'][i]))
        else:
            if (i==1):
                if (parameters.sed_flag==1):
                    mod_d.append(depth)
                    mod_d.append(depth)
                    mod_vs.append(parameters.sed_value)
                    mod_vs.append(float(modvs['vs'][i]))
                elif (parameters.sed_flag==0):
                    mod_d.append(depth)
                    mod_d.append(depth)
                    #
                    mod_vs.append(float(modvs['vs'][i-1]))
                    mod_vs.append(float(modvs['vs'][i]))
            else:
                # 2 time depth
                mod_d.append(depth)
                mod_d.append(depth)
                # 
                mod_vs.append(float(modvs['vs'][i-1]))
                mod_vs.append(float(modvs['vs'][i]))
                # the layered style of mod.STA using layer thickness. Thus start from layer 2
            thick_now = depth - allgridmoddepth[i-1]
            if (parameters.sed_flag==1):
                if (i <= sednpara) and (i > 0):
                    tmpsv0=tmpsv0+' %5.6f'%(parameters.sed_value)
                    tmpsd0=tmpsd0+' %5.6f'%((thick_now)/sedthick)
            elif (parameters.sed_flag==0):
                tmpsv0=tmpsv0+' %5.6f'%((float(modvs['vs'][i-1])+float(modvs['vs'][i]))/2)
                tmpsd0=tmpsd0+' %5.6f'%(thick_now/sedthick)
            if (i <= crustnpara+sednpara) & (i > sednpara):
                tmpcv0 += ' %5.6f'%((float(modvs['vs'][i-1]) + float(modvs['vs'][i])) / 2)
                tmpcd0=tmpcd0+' %5.6f'%(thick_now/crustthick)
            elif (i <= mantlenpara+crustnpara+2) & (i > crustnpara+sednpara+1):
                if (parameters.man_flag==1):
                    tmpmv0=tmpmv0+' %5.6f'%((float(modvs['vs'][i-1])+float(modvs['vs'][i]))/2)
                    tmpmd0=tmpmd0+' %5.6f'%(thick_now/mantlethick)
    # Set the model print line
    if (parameters.sed_flag==1):
        sthick=seddepth; 
        locals()['modsedlayer']=SedLayerId+' '+parameters.modlay+' %4.1f'%(sthick)+' %d'%(sednpara)+' %4.1f'%(float(parameters.dds))+\
        tmpsv0+' '+tmpsd0+' '+parameters.sedrhoflag+' '+parameters.sedQflag+' '+parameters.sedPflag+' '+\
        parameters.sedmodvpvs+' '+parameters.sedVpVs1+' '+parameters.sedDrho+' '+parameters.sedDrho1+' '+\
        parameters.seddVs+' '+parameters.seddvs1+' '+parameters.sedfdvs+' '+parameters.sedfdvs1
        # 
    locals()['modcrustlayer']=CrustLayerId+' '+parameters.modlay+' %4.1f %d %4.2f'%(crustthick,crustnpara,float(parameters.ddc))+tmpcv0+' '+tmpcd0+\
    ' '+parameters.crustrhoflag+' '+parameters.crustQflag+' '+parameters.crustPflag+' '+parameters.crmodvpvs+' '+\
    parameters.crVpVs1+' '+parameters.crDrho+' '+parameters.crDrho1+' '+parameters.crdVs+' '+parameters.crdvs1+\
    ' '+parameters.crfdvs+' '+parameters.crfdvs1
        # 
    if (parameters.man_flag==1):
        locals()['modmantlelayer']=ManLayerId+' '+parameters.modlay+' %4.1f %d %4.1f'%(mantlethick,mantlenpara,float(parameters.ddm))+tmpmv0+' '+tmpmd0+\
        ' '+parameters.mantlerhoflag+' '+parameters.mantleQflag+' '+parameters.mantlePflag+' '+parameters.mtmodvpvs+' '+\
        parameters.mtVpVs1+' '+parameters.mtDrho+' '+parameters.mtDrho1+' '+parameters.mtdVs+' '+parameters.mtdVs1+' '+\
        parameters.mtfdvs+' '+parameters.mtfdvs1           
# ------------------------------- >>>>>>>> write mod.{STA} here !!! <<<<<<< ----------------------------------
    fnmod= os.path.join(data_sub, f"mod.{stadata['name']}")
    print("Write the model file in layer-cake style: ",fnmod)
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
        ff3.write(str(parameters.selection_style)+"\n")
        ff3.write(str(parameters.MC_post_process_select)+"\n")
        ff3.write(str(parameters.depth_step)+"\n")
        ff3.write(str(parameters.modout_plot_type)+"\n")
        ff3.write(str(parameters.MC_QC_thres)+"\n")
    ff3.close()
 
#     # --------------------------------------------------------------------------------------------------------
#     # -----------------------------------  Figure plot for initial setup  ------------------------------------
#     # --------------------------------------------------------------------------------------------------------
    print("Plot the input model of station %s in Layer-cake type"%(stadata['name']))
    # read in Vph, H/V, and RF info from file
    if os.path.exists(gvfile):
        allgvper,mergedgv,mergedgvun=np.loadtxt(gvfile,usecols=[0,1,2],unpack=True)
    if os.path.exists(phfile):
        allphper,mergedph,mergedphun=np.loadtxt(phfile,usecols=[0,1,2],unpack=True)
    if os.path.exists(hvfile):
        hvper,hvall,hvun=np.loadtxt(hvfile,usecols=[0,1,2],unpack=True)
    if os.path.exists(rffile):
        rft,rfamp,rfunc=np.loadtxt(rffile,usecols=[0,1,2],unpack=True)
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
    # print(mod_d, mod_vs)
    # Vs
    if (parameters.sed_flag==1):
        ax1.plot(mod_vs[0:sedindex*2],mod_d[0:sedindex*2],'y*',ms=10,label='Sed Picks')
        ax1.plot(mod_vs[sedindex*2:mohoindex*2],mod_d[sedindex*2:mohoindex*2],'ko',label='crustal Picks',alpha=0.4)
    else:
        ax1.plot(mod_vs[0:mohoindex*2],mod_d[0:mohoindex*2],'ko',label='crustal Picks',alpha=0.4)
    if (parameters.man_flag==1):
        ax1.plot(mod_vs[mohoindex*2+2:-1],mod_d[mohoindex*2+2:-1],'gs',label='mantle Picks',alpha=0.4)
    ax1.plot(mod_vs,mod_d,'r--',label='Est Vs Model')
    ax1.plot(modstavs,allgridmoddepth,'k-',label='Raw Vs Model')
    
    ax1.set_ylim([0,allgridmoddepth[-1]])
    ax1.legend(loc='lower left',fontsize=12)
    
    ax1.set_ylim(ax1.get_ylim()[::-1])
    # ax1.set_xlim([0.0,5.7])
    ax1.set_ylabel('Depth (km)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax1.set_xlabel('Vs (km/s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax1.tick_params(labeltop=True)
    ax1.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
    # PHASE
    if os.path.exists(phfile):
        ax2.errorbar(allphper,mergedph,yerr=mergedphun,fmt='.:',color='r',ecolor='k',elinewidth=1.5,capthick=1.5,label='Vph',zorder=5,alpha=0.5)
    if os.path.exists(gvfile):
        ax2.errorbar(allgvper,mergedgv,yerr=mergedgvun,fmt='.:',color='b',ecolor='k',elinewidth=1.5,capthick=1.5,label='Vgv',zorder=5,alpha=0.5)
    # ax2.errorbar(phper,mergedph[0:len(phper)],yerr=mergedphun[0:len(phper)],fmt='.:',color='k',ecolor='k',elinewidth=1.5,capthick=1.5,label='ANT Vph',zorder=5,alpha=0.5)
    #ax2.errorbar(allphper,mergedph,yerr=mergedphun,fmt='.:',color='k',ecolor='r',elinewidth=1.5,capthick=1.5,label='All Vph',zorder=5,alpha=0.5)
    ax2.set_xlabel('Period (s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax2.set_ylabel('Vph (km/s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    #ax2.set_xlim([0.0,105.0])
    # ax2.set_ylim([1.5,5.0])
    # ax2.set_xlim([0,30])
    ax2.legend(loc='lower right',fontsize=12)
    ax2.tick_params(labeltop=False)
    ax2.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)

    # H/V
    if os.path.exists(hvfile):
        ax3.errorbar(hvper,hvall,yerr=hvun,fmt='.:',color='r',ecolor='k',elinewidth=1.5,capthick=1.5,label='H/V',zorder=5,alpha=0.5)
    #ax3.errorbar(hvper,hvall,yerr=hvun,fmt='.:',color='r',ecolor='k',elinewidth=1.5,capthick=1.5,label='All H/V',zorder=5,alpha=0.5)
    ax3.set_xlim([0,30])
    ax3.set_xlabel('Period (s)',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax3.set_ylabel('H/V',fontdict = {'family':'serif','color':'darkblue','size':12,'weight':'bold'})
    ax3.set_ylim([0,2.0])
    ax3.legend(loc='upper right',fontsize=12)
    ax3.tick_params(labeltop=False)
    ax3.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)

    #RF
    if os.path.exists(rffile):
        ax4.errorbar(rft,rfamp,yerr=rfunc,fmt='s:',markersize=1,color='r',ecolor='k',elinewidth=1.5,capthick=0.0,label='RFunc',zorder=5,alpha=0.5)
        # ax4.plot(rft,rfamp,'r-',zorder=7,alpha=0.96,label='RF')
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
    

    
# if __name__ == "__setup_bspline__":
#     setup()