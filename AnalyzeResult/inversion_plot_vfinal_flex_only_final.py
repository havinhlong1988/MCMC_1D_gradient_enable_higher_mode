#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on a 15:06:20 2018
Edited for Alaska on Jul 26
Updated version of inversion_plot_manymore.py


@author: u1015716
"""

# python script to automatically pull and plot all pertinent information to plot 
# plotting: posterior distribution Vs, H/V, phase dispersion
# also, plotting misfit vs iteration #

import subprocess
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import time

if len(sys.argv[:])!=3:
    print('proper usage:')
    print('python inversion_plot_vfinal.py Sta LocationOfStaDir')
    print(len(sys.argv[:]))
    sys.exit()

#%%
mdir=str(sys.argv[2])
pwd = os.getcwd()
print("current directory: ",pwd)
# scriptdir='/home/hhhuang/Research/Bayesian_joint/MCMC_AnalyzeResult'

#load and find nearest
def find_nearest(array,value):
    #may be useful to determine grid value closest to station from Hongrui Results
    idx=(np.abs(array-value)).argmin();
    return array[idx],idx

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


#%% grab things to start...


sta=str(sys.argv[1])#PIN

indatafile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/in.data_{sta}"); 
phdatafile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/{sta}.ph"); 
gvdatafile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/{sta}.gv"); 
hvdatafile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/{sta}.HV"); 
rfdatafile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/{sta}.RF"); 
modfile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/mod.{sta}"); 
connectorfile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/in.connector"); 

if os.path.exists(indatafile): 
    # print("File {} does exist! ".format(indatafile))
    data = np.genfromtxt(indatafile, dtype=int)
else:
    print("File {} does not exist! Stop!!!".format(indatafile))
    sys.exist(0)
if (len(data)==4):
    iph=data[0]; 
    igv=data[1]; 
    ihv=data[2]; 
    irf=data[3]; 
else:
    print("The {} file format wrong! 4 rows and 1 column integer".format(indatafile))
    sys.exist(0)

print("making figure for station: ",sta)
#%% Plot type for the 1D model
# plot_type = 0: layer-cake model
# plot_type = 1: gradient model (with B-spline)
# Now we can check plot_type in data/{sta}_data/in.connector
with open(connectorfile, 'r') as file:
    for i, line in enumerate(file, start=1):
        if i == 7:
            try:
                # Split the line into individual values and convert them to integers
                plot_type = [int(value) for value in line.strip().split()]
                plot_type = int(plot_type[0]); 
            except ValueError:
                # Handle the case where a value cannot be converted to an integer
                print(f"Skipping invalid value in line {i}: {line.strip()}")
            break  # Exit the loop once the target line is processed
file.close()
# print("Plot type: ",plot_type)
# time.sleep(5)
# plot_type=0 

if (plot_type==1):
    print("Plot the model in gradient style")
    # time.sleep(2)
elif (plot_type==0):
    print("Plot the model in layer-cake style")
    # time.sleep(2)
else:
    print("Error! unknown type of model plot")
    sys.exit(0)




nmodelsprint=''

nameplus=''#str(sys.argv[3])
#nameplus=''#'_round1'

nmodelsprint=file_len(mdir+'/'+sta+nameplus+'/MC.'+sta+'.all.parameters')
print("number of model: ",nmodelsprint)

posteriorfile=mdir+'/'+sta+nameplus+'/MC.'+sta+'.all.value'
posteriorfileph=mdir+'/'+sta+nameplus+'/MC.'+sta+'.all.ph'
posteriorfilerf=mdir+'/'+sta+nameplus+'/MC.'+sta+'.all.rf'
posteriorfilee=mdir+'/'+sta+nameplus+'/MC.'+sta+'.all.e'

print("posteriorfile: ",posteriorfile)

posteriorallvals=mdir+'/'+sta+nameplus+'/MC.'+sta+'.all.parameters'

modelspacemaxfile=mdir+'/'+sta+nameplus+'/MC.'+sta+'.maxmodelvalue'
modelspaceminfile=mdir+'/'+sta+nameplus+'/MC.'+sta+'.minmodelvalue'
#'''
# %% see how misfit changes over iteration number

## If ran to also look at all posterior models animated...
cmd='grep -R "test MC.C" '+mdir+'/'+sta+nameplus+'/output_'+sta+'_tmp.txt > '+mdir+'/misfitcheck_'+sta+'.txt'
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()

# print(process.returncode)

cmd='grep -R "kai_cri" '+mdir+'/'+sta+nameplus+'/output_'+sta+'.txt >> '+mdir+'/misfitcheck_'+sta+'.txt'
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
# print(process.returncode)


nsamplesvisit,misfitfn=np.genfromtxt(mdir+'/misfitcheck_'+sta+'.txt',skip_footer=1,usecols=(3,7),unpack=True)
kaimin=np.genfromtxt(mdir+'/misfitcheck_'+sta+'.txt',skip_header=len(nsamplesvisit),usecols=[1],unpack=True)
print('kaimin...: ',kaimin)
print("Starting a long-running process...")
print("It will take several minutes, do not panic!")

#%% #get format of posteriorph and posteriorHV to then read in all data
with open(posteriorfileph,'r') as ff:
    firstline=ff.readline().rstrip()
x=firstline.split()
x=[float(i) for i in x[2:]]
tmp=0
countper=0
posteriorphper=[]
for ii in np.arange(len(x)):
    if x[ii]>=tmp or x[ii]>10:
        tmp=x[ii]
        countper=countper+1
        posteriorphper.append(x[ii])
    else:
        break

posteriorph=np.loadtxt(posteriorfileph,usecols=range(countper+2,len(x)+2))

### H/V ###
with open(posteriorfilee,'r') as ff:
    firstline=ff.readline().rstrip()
x=firstline.split()
x=[float(i) for i in x[2:]]
tmp=0
countper=0
posterioreper=[]
for ii in np.arange(len(x)):
    if x[ii]>=tmp or x[ii]>10:
        tmp=x[ii]
        countper=countper+1
        posterioreper.append(x[ii])
    else:
        break

posteriore=np.loadtxt(posteriorfilee,usecols=range(countper+2,len(x)+2))



#%% read in posterior files of Receiver function
posteriorrftime=[];posteriorfval=[]
readthis=[]

startread=0;
starttmp=[]; endtmp=[]
with open(posteriorfilerf,'r') as f:
    for num, line in enumerate(f, 1):
        if ">" in line and startread==0:
            startread=1; starttmp.append(num)
            continue
            
                     
        elif ">" in line and startread==1:
            endtmp.append(num-1); #startread=0
            starttmp.append(num)
            continue
endtmp.append(num)            

for ii in np.arange(len(endtmp)):
    timetmp,rftmp=np.genfromtxt(posteriorfilerf,usecols=(0,1),skip_header=(starttmp[ii]),max_rows=(endtmp[ii])-(starttmp[ii]),unpack=True)
    posteriorrftime.append(timetmp); posteriorfval.append(rftmp)
#'''

#%% New way of getting posteriorVs
# posterior models of Vs
with open(posteriorfile,'r') as ff:
    content=ff.readlines()
posteriorVstmp=[x.strip().split() for x in content]
# posterior models of depth, correspons to the posterior models Vs
with open(posteriorfile+'_depth','r') as ff:
    content=ff.readlines()
posteriordepthtmp=[x.strip().split() for x in content]

posteriorVs=[]
posteriordepth=[]
for ii in np.arange(len(posteriorVstmp)):
    vstmp=posteriorVstmp[ii]
    vsouttmp=[]
    if (plot_type==1):
        for jj in np.arange(int(len(vstmp)/2)):
            vsouttmp.append(vstmp[jj*2])
    elif (plot_type==0):
        for jj in np.arange(int(len(vstmp))):
            vsouttmp.append(vstmp[jj])

    posteriorVs.append(np.array(vsouttmp,dtype='float'))
    
for ii in np.arange(len(posteriordepthtmp)):
    pstmp=posteriordepthtmp[ii]
    psouttmp=[]
    if (plot_type==1):
        for jj in np.arange(int(len(pstmp)/2)):
            psouttmp.append(pstmp[jj*2])
    elif (plot_type==0):
        for jj in np.arange(int(len(pstmp))):
            psouttmp.append(pstmp[jj])

    posteriordepth.append(np.array(psouttmp,dtype='float'))
    # posteriordepth = np.array(posteriordepth)
#print(posteriordepth.shape)
#'''
f = open(sta+'_Vsv.txt','w+')
for i in range(len(posteriordepth)):
	for x in range(len(posteriordepth[i])):
		f.writelines(str(posteriordepth[i][x])+' '+str(posteriorVs[i][x])+'\n')
f.close()
#f = open(sta+'_Vdepth.txt','w+')
#for i in range(len(posteriordepth)):
#	f.writelines(str(posteriordepth[i])+' '+str(posteriorVs[i])+'\n')
#f.close()
#'''
#%%
#read in min/max looked over in model space
if plot_type==0: # Layercake, take all values
    mspacemindepth,mspacemin=np.genfromtxt(modelspacemaxfile,skip_footer=6,usecols=(0,1),unpack=True)
    mspacemaxdepth,mspacemax=np.genfromtxt(modelspaceminfile,skip_footer=6,usecols=(0,1),unpack=True)
elif plot_type==1: # Gradient, take only odd values
    mspacemaxdepth,mspacemax = np.genfromtxt(modelspacemaxfile,skip_footer=6,usecols=(0,1),unpack=True)
    mspacemaxdepth = mspacemaxdepth[1::2]
    mspacemax = mspacemax[1::2]

    mspacemindepth,mspacemin = np.genfromtxt(modelspaceminfile,skip_footer=6,usecols=(0,1),unpack=True)
    mspacemindepth = mspacemindepth[1::2]
    mspacemin = mspacemin[1::2]


#%%

## read in Inversion model results
#final velocity model
avgdepth,avgvs=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.acc.average.mod',usecols=(0,1),unpack=True)

#starting velocity model
rmdepthtmp,rmvstmp,rmvptmp,rmrhotmp=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.mod.q',usecols=(0,1,2,3),unpack=True)


rmdepth=[]; rmvs=[]; rmvpt=[]; rmrhot=[]
for ii in np.arange(len(rmdepthtmp)-1):
    if (plot_type==1):
        if rmdepthtmp[ii+1] != rmdepthtmp[ii]:
            rmdepth.append(rmdepthtmp[ii])
            rmvs.append(rmvstmp[ii])
            rmvpt.append(rmvptmp[ii])
            rmrhot.append(rmrhotmp[ii])
    elif (plot_type==0):
        rmdepth.append(rmdepthtmp[ii])
        rmvs.append(rmvstmp[ii])
        rmvpt.append(rmvptmp[ii])
        rmrhot.append(rmrhotmp[ii])

rmthick=[]; rmvsf=[]; rmvp=[]; rmrho=[]
for ii in np.arange(len(rmdepth)-1):
    rmthick.append((rmdepth[ii+1]-rmdepth[ii]))
    rmvsf.append((rmvs[ii+1]+rmvs[ii])*0.5)
    rmvp.append((rmvpt[ii+1]+rmvpt[ii])*0.5)
    rmrho.append((rmrhot[ii+1]+rmrhot[ii])*0.5)
 
## PHASE-ONLY!!!!
rmvsucvm=rmvs; 
rmdepthucvm=rmdepth

avgdepthf=avgdepth
avgvsf=avgvs
    
#%% get starting model phase and receiver function

startfileph=mdir+'/'+sta+nameplus+'/Initial.ph'
startfilehv=mdir+'/'+sta+nameplus+'/Initial.hv'
startfilerf=mdir+'/'+sta+nameplus+'/Initial.rf'


### Phase Velocity Read-In ###
# with open(startfileph,'r') as ff:
#     firstline=ff.readline().rstrip()
# x=firstline.split()
# x=[float(i) for i in x]
# #x=[float(i) for i in x[2:]]
# tmp=0
# countper=0
# fcheckper=[]
# for ii in np.arange(len(x)):
#     if x[ii]>=tmp or x[ii]>10:
#         tmp=x[ii]
#         countper=countper+1
#         fcheckper.append(x[ii])
#     else:
#         break


# fcheckRC=np.loadtxt(startfileph,usecols=range(countper,len(x)))
fcheckper,fcheckRC=np.loadtxt(startfileph,usecols=(0,1),unpack='True')
### H/V Read-In ###
# with open(startfilehv,'r') as ff:
#     firstline=ff.readline().rstrip()
# x=firstline.split()
# x=[float(i) for i in x]
# #x=[float(i) for i in x[2:]]
# tmp=0
# countper=0
# fcheckperhv=[]
# for ii in np.arange(len(x)):
#     if x[ii]>=tmp or x[ii]>10:
#         tmp=x[ii]
#         countper=countper+1
#         fcheckperhv.append(x[ii])
#     else:
#         break


# fcheckHV=np.loadtxt(startfilehv,usecols=range(countper,len(x)))
fcheckperhv,fcheckHV=np.loadtxt(startfilehv,usecols=(0,1),unpack='True')
## receiver function from starting model

timerf0,rf0=np.genfromtxt(startfilerf,usecols=(0,1),unpack='True')


#
##plt.plot(timerf0,rf0)
##timerfcheck,rfcheck=np.genfromtxt('/uufs/chpc.utah.edu/common/home/u1015716/USC_Scripts/joint_inversion/MCMC_FINAL_JOINT_orig/example/data/X14A.RF',usecols=(0,1),unpack='True')
##plt.plot(timerfcheck,rfcheck)
#
##%% check forward modeling via Herrmann of MCMC code, can skip this if you want
#
#avgdepthnew,avgvsnew=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.acc.average.mod.q',usecols=(0,1),unpack=True)
#### forward model the final model average to compare and check phase velocity it outputs
#### Run forward model and save results for initial H/V, Vp, Vs for each layer ###
#def forwardmodel( velocityfile, moveddispfile ):
#    os.chdir(scriptdir)
#    # # for usage, type sprep96 -h or doc ~/PROGRAMS.330/DOC/OVERVIEW.pdf/cps330o.pdf
#    os.system('rm SREGN.ASC') #safety
#    os.system('sprep96 -M '+velocityfile+' -PARR period_input.txt -L -R') # generate preparation file for forward calculation
#    os.system('sdisp96') # compute dispersion
#    os.system('sregn96') # compute eigenfunctions
#    os.system('sdpegn96 -E -ASC -R') # this outputs in file SREGN.ASC 
#    # compute theoretical dispersion, ellipticity for Rayleigh wave
#    #  col - name as 3 - period, 4 - freq, 5 - Vph, 6 - Vgp, 9 - Ellipicitity
#    os.system('mv SREGN.ASC '+moveddispfile)
#    return;
###################################
#
#def getvp(vstmp):
#    vptmp=0.9409+(2.0947*vstmp)+(-0.8206*vstmp*vstmp)+(0.2683*vstmp*vstmp*vstmp)+(-0.0251*vstmp*vstmp*vstmp*vstmp)
#    return vptmp
#
#def getrho(vstmp):
#    rhotmp = 1.22679 + (1.53201*vstmp) -(0.83668*vstmp*vstmp) + (0.20673*vstmp*vstmp*vstmp) -(0.01656*vstmp*vstmp*vstmp*vstmp)
#    return rhotmp
##
##vp=[]; rho=[]
##for ii in np.arange(len(avgvsnew)):
##    vp.append(getvp(avgvsnew[ii])); rho.append(getrho(avgvsnew[ii]))
##    
###### Get starting model H/V and phase ####
##fcheckdepth=[]; fcheckthick=[]; fcheckvp=[]; fcheckvs=[]; fcheckrho=[]
##for ii in np.arange(len(avgdepthnew)-1):
##    index=int(ii); index1=int(ii+1)
##    if not ((avgdepthnew[index1]-avgdepthnew[index]))==0:
##        fcheckdepth.append((avgdepthnew[index1]+avgdepthnew[index])*0.5)
##        fcheckthick.append(avgdepthnew[index1]-avgdepthnew[index])
##        fcheckvp.append((vp[index]+vp[index1])*0.5)
##        fcheckvs.append((avgvsnew[index]+avgvsnew[index1])*0.5)
##        fcheckrho.append((rho[index]+rho[index1])*0.5)
##
###### get forward model for final average model to check things... ####
##Qp=np.zeros(len(fcheckthick)); Qs=Qp
##etap=np.ones(len(fcheckthick)); etas=etap; frefp=etas; frefs=etas #real
##modelfile=scriptdir+'/'+sta+'_tmp_tmodel_avg.txt'
##
### first get top 12 lines for new forward modeling code using example input_model file
##fwdscriptdir='/uufs/chpc.utah.edu/common/home/flin-group3/emberg/USC_data_2015/HV/MODEL'
##os.system('head -n12 '+fwdscriptdir+'/Herrmann/input_model > '+modelfile)
##
###insert thickness, Vp, Vs, rho, Qp=0, Qs=0, ETAP=1, ETAS=1, FREFP=1, FREFS=1
##fullprint=np.matrix([fcheckthick,fcheckvp,fcheckvs,fcheckrho,Qp,Qs,etap,etas,frefp,frefs])
##with open ( modelfile, 'a' ) as ff:
##    np.savetxt(ff, fullprint.T, fmt='%5.5f')
##
##modelfiletmp=sta+'_tmp_tmodel_avg.txt'
##forwardmodel(modelfiletmp,'orig.asc') #use this as doesn't go over allotted amount of chars in str
##periodcheck,phasecheck,hvstart=np.loadtxt('orig.asc',skiprows=1,usecols=(2,4,8),unpack=True)
##
##
#### check starting model...
##### Get starting model H/V and phase ####
##rmdepthtmp,rmvstmp
#vp=[]; rho=[]
#for ii in np.arange(len(rmvstmp)):
#    vp.append(getvp(rmvstmp[ii])); rho.append(getrho(rmvstmp[ii]))
#    
#fcheckdepth=[]; fcheckthick=[]; fcheckvp=[]; fcheckvs=[]; fcheckrho=[]
#for ii in np.arange(len(rmdepthtmp)-1):
#    index=int(ii); index1=int(ii+1)
#    if not ((rmdepthtmp[index1]-rmdepthtmp[index]))==0:
#        fcheckdepth.append((rmdepthtmp[index1]+rmdepthtmp[index])*0.5)
#        fcheckthick.append(rmdepthtmp[index1]-rmdepthtmp[index])
#        fcheckvp.append((vp[index]+vp[index1])*0.5)
#        fcheckvs.append((rmvstmp[index]+rmvstmp[index1])*0.5)
#        fcheckrho.append((rho[index]+rho[index1])*0.5)
#
#
##### get forward model of starting model to check things ####
#Qp=np.zeros(len(fcheckthick)); Qs=Qp
#etap=np.ones(len(fcheckthick)); etas=etap; frefp=etas; frefs=etas #real
#
#modelfile=scriptdir+'/tmp_fwdmodel_starting.txt'
## first get top 12 lines for new forward modeling code using example input_model file
#fwdscriptdir='/uufs/chpc.utah.edu/common/home/flin-group3/emberg/USC_data_2015/HV/MODEL'
#os.system('head -n12 '+fwdscriptdir+'/Herrmann/input_model > '+modelfile)
#
##insert thickness, Vp, Vs, rho, Qp=0, Qs=0, ETAP=1, ETAS=1, FREFP=1, FREFS=1
#fullprint=np.matrix([fcheckthick,fcheckvp,fcheckvs,fcheckrho,Qp,Qs,etap,etas,frefp,frefs])
#with open ( modelfile, 'a' ) as ff:
#    np.savetxt(ff, fullprint.T, fmt='%5.5f')
#
#forwardmodel(modelfile,'orig_start.asc')
#periodcheckstart,phasecheckstart,hvstart=np.loadtxt(scriptdir+'/orig_start.asc',skiprows=1,usecols=(2,4),unpack=True)


#%%
# final average velocity model
# if (plot_type==1):
avgdepthnew,avgvsnew,avgnewstd,avgnewstdn,avgnewstdp=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.ave.mod',usecols=(0,1,2,5,6),unpack=True)
# elif (plot_type==0):
#     avgdepthnew0,avgvsnew0,avgnewstdn,avgnewstdp=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.ave.mod',usecols=(0,1,5,6),unpack=True)
#     avgvsnew = [];
#     avgdepthnew = [];
   
#     for ivs in range(len(avgvsnew0)+0):
        # if ivs == 0:
        #     avgdepthnew.append(avgdepthnew0[ivs])
        #     avgvsnew.append(avgvsnew0[ivs])
        # else:
        #     avgdepthnew.append(avgdepthnew0[ivs])
        #     avgvsnew.append(avgvsnew0[ivs-1])
        #     avgdepthnew.append(avgdepthnew0[ivs])
        #     avgvsnew.append(avgvsnew0[ivs])

f = open(sta+'_Vsv_average.txt','w+')
print("write in!")
for i in range(len(avgdepthnew)):
	f.writelines(str(avgdepthnew[i])+' '+str(avgvsnew[i])+'\n')
f.close()
print("write in end!")

# final minmisfit velocity model
mindepthnew,minvsnew=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.minmisfit.mod',usecols=(0,1),unpack=True)

# misfit information
with open(mdir+'/'+sta+nameplus+'/MC.'+sta+'.acc.average.info','r') as misff:
    firstline=misff.readline().rstrip()
x=firstline.split()
misfittpm=np.float64(x[2])
misfitprint=str(round(misfittpm,2))

# minmisfit information
with open(mdir+'/'+sta+nameplus+'/MC.'+sta+'.minmisfit.info','r') as misff:
    firstline=misff.readline().rstrip()
x=firstline.split()
minmisfit=np.float64(x[2])
minmisfitprint=str(round(minmisfit,2))


#get depth information
depth=[]
#looping through 1 model's depth range
for ii in np.arange(0,len(posteriorVs[0])):
    depth.append(ii)
    #depth.append(ii+1)

#plot ALL models, 1std and average
#loop through all depths
avgVs=[]
stdVsincrease=[]
stdVsdecrease=[]
matrixhist=[]
minvs=2.0; maxvs=5.0

#%%

#final phase from inversion
if (iph==1):
    mcper,mcpf,mcpo,mcpounc=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.acc.average.p.disp',usecols=(0,1,2,3),unpack=True)

#final h/v from inversion
if (ihv==1):
    mceper,mchvf,mchvo,mchvounc=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.acc.average.e.disp',usecols=(0,1,2,3),unpack=True)
    mchvounc=np.loadtxt(hvdatafile,usecols=(2),unpack=True)


#final rf from inversion
mcrftime,mcrfval,mcrotime,mcrfo,mcrfounc=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.acc.average.rf',usecols=(0,1,2,3,4),unpack=True)


#minimum misfit from inversion..
mincper,mincpf,mincpo,mincpounc=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.minmisfit.p.disp',usecols=(0,1,2,3),unpack=True)
#
minceper,minchvf,minchvo,minchvounc=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.minmisfit.e.disp',usecols=(0,1,2,3),unpack=True)
#
mincrftime,mincrfval,mincrotime,mcinrfo,mincrfounc=np.loadtxt(mdir+'/'+sta+nameplus+'/MC.'+sta+'.minmisfit.rf',usecols=(0,1,2,3,4),unpack=True)



#%% Extra misfit info, put into top of Vs title plot
fname=mdir+'/'+sta+nameplus+'/MC.'+sta+'.acc.average.info'

#check if file exists (then inversion ran)
if os.path.isfile(fname):
    
    with open(fname,'r') as ff:
        firstline=ff.readline().rstrip()
    x=firstline.split()
    misfitavgmodel=float(x[2])

#%% Get thickness variations
# if there is a min/max parameters file, get min/max sed and crust thickness and plot
# test for sed/crust thickness change combination first and then reapply as necessary
# original thickness: sta+'_data/mod.ucvm.'+sta  , 

# with open(modfile) as fp:
#     for i, line in enumerate(fp):
#         if i == 0:
#             x=line.split()
#             sthick0=float(x[2])
#         if i == 1:
#             x=line.split()
#             cthick0=sthick0+float(x[2])
# cthickmin=0
# fname=mdir+'/'+sta+nameplus+'/MC.'+sta+'.minparas'            
# with open(fname,'r') as ff:
#     firstline=ff.readline().rstrip()
# x=firstline.split()
# sthickmin=float(x[len(x)-2]) #comment if not changing crust thick
# cthickmin=float(x[len(x)-1])

# #%%
# cthickmax=0
# fname=mdir+'/'+sta+nameplus+'/MC.'+sta+'.maxparas'            
# with open(fname,'r') as ff:
#     firstline=ff.readline().rstrip()
# x=firstline.split()
# sthickmax=float(x[len(x)-2]) #comment if not changing crust thick
# cthickmax=float(x[len(x)-1])
#%%
sthickf=0
fname=mdir+'/'+sta+nameplus+'/MC.'+sta+'.ave.paras'            
with open(fname,'r') as ff:
    firstline=ff.readline().rstrip()
x=firstline.split()
sthickf=float(x[len(x)-2]) #comment if not changing crust thick
cthickf=sthickf+float(x[len(x)-1])


#MC.A19K.minmisfit.mod.group , formatted like model input file...
with open(mdir+'/'+sta+nameplus+'/MC.'+sta+'.minmisfit.mod.group','r') as ff:
    datatmp=ff.readlines()
# crustminf=float(datatmp[0].split()[2])+float(datatmp[1].split()[2])
# sthickminf=float(datatmp[0].split()[2])

#%% New plotting AND can pick coordinates to run over for individual points, OR for cross-section :)
#from mpl_toolkits.basemap import Basemap, addcyclic

print('making figure...')

#fig=plt.figure(1,figsize=(13,18))
#
#ax1=plt.subplot2grid((3,2),(1,0),rowspan=2)
##ax1=plt.subplot2grid((3,2), (0,0), rowspan=2)
#ax2=plt.subplot2grid((3,2), (0,1))
#ax3=plt.subplot2grid((3,2), (1,1))
#ax4=plt.subplot2grid((3,2), (0,0)) #misfit
#
##ax4=plt.subplot2grid((3,2), (2,0)) #misfit
#ax5=plt.subplot2grid((3,2), (2,1))

### NEW!!
fig=plt.figure(int(nmodelsprint),figsize=(18,20))

#Alt 1
#ax1=plt.subplot2grid((3,3),(1,0),rowspan=2)
##ax1=plt.subplot2grid((3,2), (0,0), rowspan=2)
#ax2=plt.subplot2grid((3,3), (0,1),colspan=2)
#ax3=plt.subplot2grid((3,3), (1,1),colspan=2)
#ax4=plt.subplot2grid((3,3), (0,0)) #misfit
#ax5=plt.subplot2grid((3,3), (2,1),colspan=2)

#Alt 2
ax1=plt.subplot2grid((3,5),(1,0),rowspan=2,colspan=2)
ax2=plt.subplot2grid((3,5), (0,2),colspan=3)
ax3=plt.subplot2grid((3,5), (1,2),colspan=3)
ax4=plt.subplot2grid((3,5), (0,0),colspan=2)
ax5=plt.subplot2grid((3,5), (2,2),colspan=3)

ax1.set_facecolor('lightyellow')
ax2.set_facecolor('lightyellow')
ax3.set_facecolor('lightyellow')
ax4.set_facecolor('lightyellow')
ax5.set_facecolor('lightyellow')


#Vs


#plot all models via colors
#color=iter(plt.cm.rainbow(np.linspace(0,1,len(posteriorVs))))
for jj in np.arange(len(posteriorVs)-1):
    #c=next(color)
    c='lightgray'
    ax1.plot(posteriorVs[jj],posteriordepth[jj],c=c,alpha=0.1)
    
# ax1.plot(minvsnew,mindepthnew,'bs--',lw=2,ms=10,label='MinMisfit Vs',zorder=7,mec='k') # Min misfit model
ax1.plot(avgvsnew,avgdepthnew,'r',linewidth=3,ms=10,label='Final Vs',zorder=4,alpha=0.5) # Average model
# ax1.errorbar(avgvsnew,avgdepthnew, xerr=avgnewstd, fmt='-o',ecolor="m",elinewidth=2,capsize=4,alpha=1.0)
# ax1.plot(avgvsnew,avgdepthnew,'yo',lw=0.5,ms=5,label='Mean Vs',zorder=7,mec='k')
#add back in after including mantle area searched!!!!
# ax1.plot(mspacemin,mspacemindepth,'g',label='ModelSpace',dashes=[8,8,16,8]) #on, off, on, off
# ax1.plot(mspacemax,mspacemaxdepth,'g',dashes=[8,8,16,8])
#ax1.plot(avgnewstdp,avgdepthnew,'k--',label='StdDev')
#ax1.plot(avgnewstdn,avgdepthnew,'k--')

## thickness stuff
# ax1.plot([0,8],[sthickmin,sthickmin],'k--',lw=1.5,zorder=1,label='sed min_dep') #comment if not changing crust thick
# ax1.plot([0,8],[sthickmax,sthickmax],'b--',lw=1.5,zorder=1,label='sed max_dep') #commint if not changing crust thick
##
#ax1.plot([0,5],[cthickmin,cthickmin],'b--',lw=0.5,zorder=1)
#ax1.plot([0,5],[cthickmax,cthickmax],'b--',lw=0.5,zorder=1)
##
#ax1.plot([0,5],[sthick0,sthick0],'k.-',alpha=0.5,zorder=1)
#ax1.plot([0,5],[cthick0,cthick0],'k.-',alpha=0.5,zorder=1)
#
#ax1.plot([0,5],[cthickf,cthickf],'c-.',lw=3,zorder=2)
#ax1.plot([0,5],[sthickf,sthickf],'c-.',lw=3,zorder=2) #comment if not changing crust thick
#
#ax1.plot([0,5],[crustminf,crustminf],'g-.',lw=3,zorder=2)
#ax1.plot([0,5],[sthickminf,sthickminf],'g-.',lw=3,zorder=2)


###
# ax1.plot(rmvsucvm,rmdepthucvm,'r^-',lw=2.5,ms=10,label='Start Vs')

ax1.legend(loc='lower left',fontsize=15)
#ax1.set_yticks(np.arange(0, int(round(np.max(rmdepthucvm),0)), 25.0))
ax1.set_yticks(np.arange(0, 50.0, 5.0))
#ax12=plt.gca()
ax1.set_xlabel('Vs (km/s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax1.set_ylabel('Depth (km)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.tick_params(labeltop=False,labelsize=20)#, labelright=True)
#ax1.set_ylim([0.0, 150])
ax1.set_ylim([0.0, 60])
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_xlim([0.0,7.0])
ax1.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)

print('Plot inversion output')
title='Vs of station: ['+sta +'] - MinimumMisfit: '+str(round(minmisfit,2))+' - #posteriorModels:'+str(len(posteriorfval))
ax4.set_title(title,y=1.1,loc='left',fontdict = {'family':'serif','color':'darkgreen','size':30,'weight':'bold'}) #used to be 1.04

# Vs Upper 20(?) km
#title='Vs   Misfit: '+str(round(misfitavgmodel,2))
#ax4.set_title(title,y=1.06) #used to be 1.04

#add back in after including mantle area searched!!!!
# ax4.plot(mspacemin,mspacemindepth,'g',label='ModelSpace',dashes=[8,8,16,8]) #on, off, on, off
# ax4.plot(mspacemax,mspacemaxdepth,'g',dashes=[8,8,16,8])
# ax4.plot([0,8],[sthickmin,sthickmin],'k--',lw=1.5,zorder=1,label='sed min_dep') #comment if not changing crust thick
# ax4.plot([0,8],[sthickmax,sthickmax],'b--',lw=1.5,zorder=1,label='sed max_dep') #commint if not changing crust thick

for jj in np.arange(len(posteriorVs)-1):
    #c=next(color)
    c='lightgray'
    ax4.plot(posteriorVs[jj],posteriordepth[jj],c=c,alpha=0.1)

#ax4.plot(posteriorVs[minmisfitpostidx],posteriordepth[minmisfitpostidx],'mo',ms=10,lw=2,label='MinMisfit Vs',zorder=7) #NEW

# ax4.plot(minvsnew,mindepthnew,'bs--',lw=2,ms=10,label='MinMisfit Vs',zorder=7,mec='k') # Min misfit model
ax4.plot(avgvsnew,avgdepthnew,'r',linewidth=3,ms=10,label='Final Vs',zorder=4,alpha=0.5)
# ax4.errorbar(avgvsnew,avgdepthnew, xerr=avgnewstd, fmt='-o',ecolor="m",elinewidth=2,capsize=4,alpha=1.0)

# ax4.plot(avgvsnew,avgdepthnew,'yo',lw=1,ms=5,label='Mean Vs',zorder=7,mec='k')

# ax4.plot(posteriorVs[minmisfitpostidx],posteriordepth[minmisfitpostidx],'yo',lw=2,label='MinMisfit Vs',zorder=7)
#ax4.plot(avgnewstdp,avgdepthnew,'k--',label='StdDev')
#ax4.plot(avgnewstdn,avgdepthnew,'k--')
# ax4.plot(rmvsucvm,rmdepthucvm,'r^-',lw=2.5,ms=10,label='Start Vs')
ax4.set_yticks(np.arange(0, int(round(np.max(rmdepthucvm),0)), 5.0))
#ax12=plt.gca()
ax4.set_xlabel('Vs (km/s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax4.set_ylabel('Depth (km)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax4.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax4.tick_params(labeltop=True,labelsize=20)#, labelright=True)
ax4.set_ylim([0.0, 5])
ax4.set_ylim(ax4.get_ylim()[::-1])
ax4.set_xlim([0.0,7.0])
ax4.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
# ax4.legend(loc='lower left',fontsize=15)

# Vp_Phase
ax3.set_title('Phase Velocity', loc='left',fontdict = {'family':'serif','color':'black','size':25,'weight':'bold'})# Misfit:'+misfitprint)
if (iph==1):
    # #color=iter(plt.cm.rainbow(np.linspace(0,1,len(posteriorph[:,0]))))
    if np.all(posteriorph[:, 0] == 0):
        print("No posterior of Vph???")
    else:
        for ii in range(len(posteriorph[:, 0]) - 1):
            c = 'lightgray'
            ax3.plot(posteriorphper, posteriorph[ii], c=c, alpha=0.5)


    ax3.errorbar(mcper,mcpo,yerr=mcpounc,fmt='o',color='k',ecolor='k',ms=7,elinewidth=2.5,capthick=1.5,label='Data Vph',zorder=6)
    # ax3.plot(mincper,mincpf,'bs',lw=1,ms=6,label='Minmisfit Vph',zorder=7,mec='k')
    ax3.plot(mcper,mcpf,'r',linewidth=3,ms=10,label='Final Vph',zorder=4,alpha=0.5)
    # ax3.plot(mcper,mcpf,'wo',lw=1,ms=6,label='Final Vph',zorder=7,mec='k')
    #ax3.plot(mcper,mcpf,'yo',lw=0.5,ms=4,label='MC Vph',zorder=7,mec='k')

# ax3.plot(fcheckper,fcheckRC,'r^',label='Starting Vph',ms=10,zorder=5)

#check forward modeling results..
#ax3.plot(periodcheck,phasecheck,'m.',zorder=7) #check final
#ax3.plot(periodcheckstart,phasecheckstart,'g.',zorder=7) #check start

#ax3.legend(loc='upper left',fontsize=15)
ax3.legend(loc='best',fontsize=15)

ax3.set_xlabel('Period (s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
#ax3.set_xlim([5.5,16.5])
ax3.set_ylabel('Vph (km/s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
#ax3.set_xlim([8,40])
# ax3.set_xlim([2.5,10.5])
# ax3.set_ylim(1.95, 4.6)#(1.95, 3.1)#
ax3.set_ylim(1.5, 5)#(1.95, 3.1)#
ax3.set_xlim([posteriorphper[0]-0.5,posteriorphper[-1]+0.5])
ax3.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
ax3.tick_params(labelsize=20)#, labelright=True)
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")

# H/V
ax5.set_title('H/V', loc='left', fontdict = {'family':'serif','color':'black','size':25,'weight':'bold'})# Misfit:'+misfitprint)
color=iter(plt.cm.rainbow(np.linspace(0,1,len(posteriore[:,0]))))
if (ihv==1):
    for ii in np.arange(len(posteriore[:,0])-1):
        c='lightgray'
        ax5.plot(posterioreper,posteriore[ii],c=c,alpha=0.5)
    ax5.errorbar(mceper,mchvo,yerr=mchvounc,fmt='o',color='k',ms=7,label='Data H/V',ecolor='k',elinewidth=2.5,capthick=1.5,zorder=6)

    ##Joint inversion plotting
    # ax5.plot(minceper,minchvf,'bs',lw=1, ms=7,label='minmisfit H/V',zorder=7,mec='k')
    # ax5.plot(mceper,mchvf,'wo',lw=1, ms=7,label='Final H/V',zorder=7,mec='k')
    ax5.plot(mceper,mchvf,'r',linewidth=3,ms=10,label='Final H/V',zorder=4,alpha=0.5)
    #ax5.plot(mceper,mchvf,'yo',lw=0.5, ms=4,label='mean MC H/V',zorder=7,mec='k')
    #ax5.plot(periodcheckstart,hvstart,'m.',zorder=7)##hvstart

# ax5.plot(fcheckperhv[:len(fcheckHV)],fcheckHV,'r^',ms=10,label='Starting H/V',zorder=5)
#ax5.plot(fcheckperhv,fcheckHV,'r^',ms=10,label='Starting H/V',zorder=5)
ax5.legend(loc='upper right',fontsize=15)

ax5.set_xlabel('Period (s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax5.set_ylabel('H/V',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax5.set_ylim([0.5, 1.5])
ax5.tick_params(labelsize=20)#, labelright=True)
ax5.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
ax5.yaxis.tick_right()
ax5.yaxis.set_label_position("right")



##############
# RF
ax2.set_title('Receiver function', loc='left',fontdict = {'family':'serif','color':'black','size':25,'weight':'bold'})# Misfit'+misfitprint)
if (irf==1):
    ax2.errorbar(mcrotime,mcrfo,yerr=mcrfounc,fmt='.:',color='k',ecolor='k',elinewidth=1.5,capthick=0.0,label='Data RFs',zorder=5,alpha=0.5)


    # ax2.plot(mincrftime,mincrfval,'bs',lw=1, ms=7,mec='k',label='Minmisfit RF',zorder=7)

    # ax2.plot(mcrftime,mcrfval,'wo',lw=1,ms=7,label='Final RF',zorder=7,mec='k')
    ax2.plot(mcrftime,mcrfval,'r',linewidth=3,ms=10,label='Final RF',zorder=4,alpha=0.5)
    # ax2.plot(mcrftime,mcrfval,'ko',ms=7,zorder=7)


    # ax2.plot(mcrftime,mcrfval,'k.',ms=6,zorder=6)
    #color=iter(plt.cm.rainbow(np.linspace(0,1,len(posteriorfval))))
    for ii in np.arange(len(posteriorfval)-1):
        c='lightgray'
        ax2.plot(posteriorrftime[ii],posteriorfval[ii],c=c,alpha=0.2)


# ax2.plot(timerf0,rf0,'r',linewidth=3,ms=10,label='Starting RF',zorder=4,alpha=0.5)
# ax2.plot(timerf0,rf0,'r^',ms=10,label='Starting RFs',zorder=5)
ax2.legend(loc='upper right',fontsize=15)

#ax2.set_ylim([min(rf0)-0.1,max(rf0)+0.1])
ax2.set_ylim([-0.25,0.6])
ax2.set_xlabel('Time (s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax2.set_ylabel('RF Amp',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax2.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
ax2.tick_params(labelsize=20)#, labelright=True)
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")

plt.savefig(mdir+'/'+sta+nameplus+'/'+sta+nameplus+'_MCMC.png',bbox_inches='tight',transparent=False,pad_inches=0.1, dpi=300)
#plt.savefig(mdir+'/'+sta+nameplus+'/'+sta+nameplus+'_MCMC.pdf',bbox_inches='tight',transparent=True,pad_inches=0.1, dpi=300,format='pdf')
plt.close()
#plt.show()

# print('Plot misfit toward iterations...')
# #%%
# fig=plt.figure(2,figsize=(5,5))

# ax1=plt.subplot(211)
# ax2=plt.subplot(212)


# #misfit function top part
# ax1.set_title('Misfit Functional: '+sta+' #misfitcrit='+str(kaimin)+' #posteriorModels:'+str(len(posteriorfval)),fontdict = {'family':'serif','color':'darkgreen','size':11,'weight':'bold'})
# ax1.axhline(kaimin,color='k',linestyle='--',lw=2)
# for ii in np.arange(len(misfitfn)):
#     if misfitfn[ii]<=kaimin:
#         ax1.plot(nsamplesvisit[ii],misfitfn[ii],'ko',alpha=0.2,ms=5,mec='k')
#     else:
#         ax1.plot(nsamplesvisit[ii],misfitfn[ii],'ro',alpha=0.2,ms=4,mec='k')
        

# #ax1.set_xlabel('Iteration #')
# ax1.set_ylabel('Zoom Misfit',fontdict = {'family':'serif','color':'darkblue','size':11,'weight':'bold'})
# ax1.set_ylim([kaimin*2,0])




# ## Misfit functional
# ax2.axhline(kaimin,color='k',linestyle='--',lw=2,label='critical misfit')
# for ii in np.arange(len(misfitfn)):
#     if misfitfn[ii]<=kaimin:
#         ax2.plot(nsamplesvisit[ii],misfitfn[ii],'ko',alpha=0.2,ms=5,mec='k')
#     else:
#         ax2.plot(nsamplesvisit[ii],misfitfn[ii],'ro',alpha=0.2,ms=4,mec='k')
        

# ax2.set_xlabel('Iteration number',fontdict = {'family':'serif','color':'darkblue','size':11,'weight':'bold'})
# ax2.set_ylabel('Misfit',fontdict = {'family':'serif','color':'darkblue','size':11,'weight':'bold'})
# ax2.set_ylim(ax2.get_ylim()[::-1])

# plt.savefig(mdir+'/'+sta+nameplus+'/'+sta+nameplus+'_IterVsMisfit.png',bbox_inches='tight',transparent=False,pad_inches=0.1, dpi=300)

# plt.close()
# #plt.show()