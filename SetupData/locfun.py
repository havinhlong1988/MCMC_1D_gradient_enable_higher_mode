import numpy as np
import os, sys
import matplotlib.pyplot as plt
import subprocess
import time 

#function to linearly interpolate
def lininterpwiki(x,x0,y0,x1,y1):
    tmp=((x-x0)/(x1-x0))
    y=(y0*(1-tmp))+(y1*tmp)
    return y

#load and find nearest
def find_nearest(array,value):
    #may be useful to determine grid value closest to station from Hongrui Results
    idx=(np.abs(array-value)).argmin();
    return array[idx],idx
#Example: lonsta,ii=find_nearest(llons,lontmp)

## Peak Detection ##
def peakdet2(v, depth):
    """
    Elizabeth Berg's version of detecting peak prominence
    Should work better (maybe?) than peakdet from Eli Billauer
    
    Similar to FTAN, if a peak is 'prominent', then include
    
    """
    if len(v)!=len(depth):
        #print len(v)
        #print '!='
        #print len(depth)
        print('ERROR!!!! len(v)!=len(depth)')
        sys.exit()
        
    allprom=[]; depout=[]
    for ii in np.arange(len(v)-2):
        h1=depth[ii+1]-depth[ii]
        h2=depth[ii+2]-depth[ii+1]
        h3=h1+h2
        
        tmp=((v[ii]/h1)-(((1./h1)+(1./h2))*v[ii+1])+(v[ii+2]/h2))*(h3/4.)*(10.)
        allprom.append(tmp)
        depout.append(depth[ii+1])
        
    return allprom,depout
