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
import pandas as pd
import time

if len(sys.argv[:])!=3:
    print('proper usage:')
    print('python inversion_plot_vfinal.py Sta LocationOfStaDir')
    print(len(sys.argv[:]))
    sys.exit()

#%%
mdir=str(sys.argv[2])
pwd = os.getcwd()
# mdir=os.path.join(pwd)
print("current directory: ",pwd)
# scriptdir='/home/hhhuang/Research/Bayesian_joint/MCMC_AnalyzeResult'
pardir = os.path.dirname(pwd)
velmodfile = os.path.join(pardir,"Vel_mod","NVN_1D_modified.dat")
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

# sta="00660"
sta=str(sys.argv[1])#PIN

indatafile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/in.data_{sta}"); 
phdatafile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/{sta}.ph"); 
gvdatafile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/{sta}.gv"); 
hvdatafile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/{sta}.HV"); 
rfdatafile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/{sta}.RF"); 
modfile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/mod.{sta}"); 
connectorfile=os.path.join(os.path.join(os.path.dirname(pwd),"data"),f"{sta}_data/in.connector"); 
#
phdataprintfile=os.path.join(pwd,sta,"Initial.ph");
hvdataprintfile=os.path.join(pwd,sta,"Initial.hv");
rfdataprintfile=os.path.join(pwd,sta,"Initial.rf");

if os.path.exists(indatafile): 
    # print("File {} does exist! ".format(indatafile))
    data = np.genfromtxt(indatafile, dtype=int)
else:
    print("File {} does not exist! Stop!!!".format(indatafile))
    sys.exist(0)
if (len(data)==9):
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
        if i == 9:
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

# ========================================= Initial Model setting =====================================================
inputmodelfile = mdir+'/'+sta+nameplus+'/Initial.mod'
inputmodelfileleft = mdir+'/'+sta+nameplus+'/Initial_left.mod'
inputmodelfileright = mdir+'/'+sta+nameplus+'/Initial_right.mod'

inputmodel0 = pd.read_csv(velmodfile,sep="\s+",usecols=[0,1,7],names=["dep","vs","vp"],header=None)
inputmodel = pd.read_csv(inputmodelfile,sep="\s+",usecols=[0,1,7],names=["dep","vs","id"],header=None)
inputmodelleft = pd.read_csv(inputmodelfileleft,sep="\s+",usecols=[0,1,7],names=["dep","vs","id"],header=None)
inputmodelright = pd.read_csv(inputmodelfileright,sep="\s+",usecols=[0,1,7],names=["dep","vs","id"],header=None)
#
# Duplicate row 0
first_row = inputmodel.iloc[0].copy()
# Change depth to 0.0
first_row["dep"] = 0.0
# Insert to the top
inputmodel = pd.concat([pd.DataFrame([first_row]), inputmodel], ignore_index=True)
# -----------------
first_row = inputmodelleft.iloc[0].copy()
# Change depth to 0.0
first_row["dep"] = 0.0
# Insert to the top
inputmodelleft = pd.concat([pd.DataFrame([first_row]), inputmodelleft], ignore_index=True)
# -----------------
first_row = inputmodelright.iloc[0].copy()
# Change depth to 0.0
first_row["dep"] = 0.0
# Insert to the top
inputmodelright = pd.concat([pd.DataFrame([first_row]), inputmodelright], ignore_index=True)
# -------------------- find moho depth ------------------------------------------------------
# Compute velocity difference between layers
inputmodel["dvs"] = inputmodel["vs"].diff()

# Find index of maximum positive jump
moho_idx = inputmodel["dvs"].idxmax()-1 # right above the jump

# Output Moho depth
moho_depth = inputmodel.loc[moho_idx, "dep"]

print("Moho depth =", moho_depth)
# ================================================== Starting forward model data =====================================================

startfileph=mdir+'/'+sta+nameplus+'/Initial.ph'
startfilehv=mdir+'/'+sta+nameplus+'/Initial.hv'
startfilerf=mdir+'/'+sta+nameplus+'/Initial.rf'

# ------------------------------------------
inputph = pd.read_csv(startfileph,sep="\s+",usecols=[0,1],names=["per","vs"],header=None)
inputhv = pd.read_csv(startfilehv,sep="\s+",usecols=[0,1],names=["per","vs"],header=None)
inputrf = pd.read_csv(startfilerf,sep="\s+",usecols=[0,1],names=["time","amp"],header=None)
# ------------------------------------------
dataph = pd.read_csv(phdatafile,sep="\s+",usecols=[0,1,2],names=["per","vs","unc"],header=None)
datahv = pd.read_csv(hvdatafile,sep="\s+",usecols=[0,1,2],names=["per","vs","unc"],header=None)
datarf = pd.read_csv(rfdatafile,sep="\s+",usecols=[0,1,2],names=["time","amp","unc"],header=None)

#'''
# %% see how misfit changes over iteration number

print('making figure...')

### NEW!!
fig=plt.figure(figsize=(18,20))

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

ax1.plot(inputmodelleft['vs'],inputmodelleft['dep'],'g',label='ModelSpace',dashes=[8,8,16,8])
ax1.plot(inputmodelright['vs'],inputmodelright['dep'],'g',dashes=[8,8,16,8])
ax1.plot(inputmodel['vs'],inputmodel['dep'],'r^-',lw=2.5,ms=10,label='Start Vs')
ax1.plot(inputmodel0['vs'],inputmodel0['dep'],'bo-',lw=2.5,ms=10,alpha=0.5,label='Input model')
ax1.hlines(moho_depth, xmin=0, xmax=8, colors='k', linestyles='--', linewidth=2.5, label='Moho depth')
ax1.text(1.5, moho_depth*1.05, f"{moho_depth:.1f} km", va='center', fontsize=33)
#
ax1.legend(loc='lower left',fontsize=15)
#ax1.set_yticks(np.arange(0, int(round(np.max(rmdepthucvm),0)), 25.0))
ax1.set_yticks(np.arange(0, 50.0, 5.0))
#ax12=plt.gca()
ax1.set_xlabel('Vs (km/s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax1.set_ylabel('Depth (km)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.tick_params(labeltop=False,labelsize=20)#, labelright=True)
#ax1.set_ylim([0.0, 150])
ax1.set_ylim([0.0, 50])
ax1.set_ylim(ax1.get_ylim()[::-1])
ax1.set_xlim([0.0,7.0])
ax1.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)

print('Plot inversion output')
title='Data of station: {}'.format(sta)
ax4.set_title(title,y=1.1,loc='left',fontdict = {'family':'serif','color':'darkgreen','size':30,'weight':'bold'}) #used to be 1.04

# Vs Upper 20(?) km
ax4.plot(inputmodelleft['vs'],inputmodelleft['dep'],'g',label='ModelSpace',dashes=[8,8,16,8])
ax4.plot(inputmodelright['vs'],inputmodelright['dep'],'g',dashes=[8,8,16,8])
ax4.plot(inputmodel['vs'],inputmodel['dep'],'r^-',lw=2.5,ms=10,label='Start Vs')
ax4.plot(inputmodel0['vs'],inputmodel0['dep'],'bo-',lw=2.5,ms=10,alpha=0.5,label='Input model')

# ax4.set_yticks(np.arange(0, int(round(np.max(rmdepthucvm),0)), 5.0))
#ax12=plt.gca()
ax4.set_xlabel('Vs (km/s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax4.set_ylabel('Depth (km)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax4.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax4.tick_params(labeltop=True,labelsize=20)#, labelright=True)
ax4.set_ylim([0.0, 5])
ax4.set_yticks(np.arange(0, 5.0, 1.0))
ax4.set_ylim(ax4.get_ylim()[::-1])
ax4.set_xlim([0.0,7.0])
ax4.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
# ax4.legend(loc='lower left',fontsize=15)

# Vp_Phase
ax3.set_title('Phase Velocity', loc='left',fontdict = {'family':'serif','color':'black','size':25,'weight':'bold'})# Misfit:'+misfitprint)

ax3.plot(inputph['per'],inputph['vs'],'r^',label='Starting Vph',ms=10,zorder=5)
ax3.errorbar(dataph['per'],dataph['vs'], yerr=dataph['unc'],fmt='o',color='k',ecolor='k',ms=7,elinewidth=1.5,capthick=1.5,label='Data Vph',zorder=6)
ax3.legend(loc='best',fontsize=15)

ax3.set_xlabel('Period (s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax3.set_ylabel('Vph (km/s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax3.set_ylim(0, 5)#(1.95, 3.1)#
ax3.set_xlim([inputph['per'].min()-0.5,inputph['per'].max()+0.5])
ax3.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
ax3.tick_params(labelsize=20)#, labelright=True)
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")

# H/V
ax5.set_title('H/V', loc='left', fontdict = {'family':'serif','color':'black','size':25,'weight':'bold'})# Misfit:'+misfitprint)
ax5.plot(inputhv['per'],inputhv['vs'],'r^',ms=10,label='Starting H/V',zorder=5)
ax5.errorbar(datahv['per'],datahv['vs'], yerr=datahv['unc'],fmt='o',color='k',ecolor='k',ms=7,elinewidth=1.5,capthick=1.5,label='Data H/V',zorder=6)
ax5.legend(loc='upper right',fontsize=15)

ax5.set_xlabel('Period (s)',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax5.set_ylabel('H/V',fontdict = {'family':'serif','color':'darkblue','size':20,'weight':'bold'})
ax5.set_ylim([0.0, 2.0])
ax5.tick_params(labelsize=20)#, labelright=True)
ax5.grid(color='k', axis='both',linestyle='-', linewidth=2,alpha=0.1,zorder=2)
ax5.yaxis.tick_right()
ax5.yaxis.set_label_position("right")

##############
# RF
ax2.set_title('Receiver function', loc='left',fontdict = {'family':'serif','color':'black','size':25,'weight':'bold'})# Misfit'+misfitprint)

# ax2.plot(timerf0,rf0,'r',linewidth=3,ms=10,label='Starting RF',zorder=4,alpha=0.5)
ax2.plot(inputrf['time'],inputrf['amp'],'r-',ms=10,label='Starting RFs',zorder=5)
ax2.errorbar(datarf['time'],datarf['amp'], yerr=datarf['unc'],fmt='o',color='k',ecolor='k',ms=3,elinewidth=1.0,capthick=0.5,label='Data RF',zorder=6)
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