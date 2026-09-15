#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-project overrides for the FM project (loaded automatically by
SetupData/main.py when it is run as:  python main.py FM).

Only the names below differ from SetupData/parameters.py; everything else is
inherited.  CHT is unaffected -- it has no parameters_CHT.py, so it keeps using
parameters.py exactly as before.
"""

# --- starting model -----------------------------------------------------
# FM data are crustal-scale (Vph 3-40 s, RF to 7.5 s) and need the 50 km
# starting model, not the 15 km one CHT uses.
modfile = 'NTW1d_H14_intp_50km'

# Crust sub-layer thickness.  The FM crust is ~30 km thick, so 0.5 km
# sub-layers give ~60 layers; CHT's 0.05 km would give 600 and blow past the
# forward model's sub-layer limit.
ddc = 0.5

# --- data weights -------------------------------------------------------
# FM inverts Vph + H/V + RF jointly.  parameters.py ships phw=1.0 with the
# others at 0, which would read and plot HV/RF but give them ZERO weight in
# the misfit.  These are the weights the reference project used.
phw = 0.335
gvw = 0.0
hvw = 0.33
rfw = 0.335

# --- receiver function --------------------------------------------------
# MUST match the gaussian the observed RFs were computed with.  The FM RFs
# were made with 3.0; CHT's default of 2.5 would bias the RF fit.
RF_gaussian_width = 3.0

# --- run size -----------------------------------------------------------
# 1 thread: the Fortran forward model is called without the OpenMP critical
# section (MC.C), so multi-threaded runs can corrupt or crash on macOS.
# Raise this on CentOS, or once that race is fixed.
MC_number_of_cores = 1
