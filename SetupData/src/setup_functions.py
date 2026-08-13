#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Imports...
import numpy as np
import os, sys
import matplotlib.pyplot as plt
import subprocess
import time 
############################################## function to use in this script ##################################################### 
#function to linearly interpolate
def lininterpwiki(x,x0,y0,x1,y1):
    tmp=((x-x0)/(x1-x0))
    y=(y0*(1-tmp))+(y1*tmp)
    return y
# ----------------------------------------------------------------------------------
#load and find nearest
def find_nearest(array,value):
    #may be useful to determine grid value closest to station from Hongrui Results
    idx=(np.abs(array-value)).argmin();
    return array[idx],idx
#Example: lonsta,ii=find_nearest(llons,lontmp)
# ----------------------------------------------------------------------------------
##%% take derivative, find peaks, create starting model...
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

def detect_moho_from_layer_cake(dep, vs, 
                                min_depth=5.0,   # km, ignore very shallow jumps
                                max_depth=80.0,  # km, ignore too deep
                                min_jump=0.2):   # km/s, minimum Vs jump to be considered Moho
    """
    Detect Moho from a 1-D layer-cake model using Vs jump.

    Parameters
    ----------
    dep : 1D array-like
        Depth array (km), positive downward.
    vs : 1D array-like
        Shear-wave velocity array (km/s), same length as dep.
    min_depth, max_depth : float
        Only consider jumps within this depth range.
    min_jump : float
        Minimum Vs jump (km/s) to accept as Moho.

    Returns
    -------
    moho_depth : float
        Estimated Moho depth (km).
    idx : int
        Index i such that Moho lies between dep[i] and dep[i+1].
    """

    dep = np.asarray(dep)
    vs  = np.asarray(vs)

    # Sanity check: ensure monotonic increasing depth
    if not np.all(np.diff(dep) > 0):
        raise ValueError("Depth array must be strictly increasing.")

    # Compute layer-to-layer Vs jump
    d_vs = vs[1:] - vs[:-1]

    # Depth at top of each jump
    dep_top = dep[:-1]

    # Only consider jumps in the desired depth window
    mask = (dep_top >= min_depth) & (dep_top <= max_depth)

    # Set jumps outside the window to -inf so they are never picked
    d_vs_masked = np.where(mask, d_vs, -np.inf)

    # Find index of the largest positive Vs jump in the window
    idx = int(np.argmax(d_vs_masked))
    max_jump_val = d_vs_masked[idx]

    if max_jump_val < min_jump:
        raise RuntimeError(
            f"No Vs jump larger than {min_jump:.2f} km/s found between "
            f"{min_depth} and {max_depth} km."
        )

    # Put Moho depth at the middle of the jump (between two layers)
    moho_depth = 0.5 * (dep[idx] + dep[idx + 1])

    return moho_depth, idx

def detect_moho_duplicate_depth(depth, vs, min_jump=0.2):
    """
    Detect Moho when it is encoded as two consecutive points
    with the same depth but different Vs.

    Parameters
    ----------
    depth : 1D array-like
        Depth array (km).
    vs : 1D array-like
        Vs array (km/s).
    min_jump : float
        Minimum Vs jump (km/s) to consider as Moho.

    Returns
    -------
    moho_depth : float
    idx : int
        Index i such that Moho is between (i, i+1).
    """
    depth = np.asarray(depth)
    vs = np.asarray(vs)

    for i in range(len(depth) - 1):
        if np.isclose(depth[i], depth[i+1]):
            dv = vs[i+1] - vs[i]
            if dv > min_jump:  # positive jump big enough
                return depth[i], i

    raise RuntimeError("No Moho-like duplicate-depth jump found.")