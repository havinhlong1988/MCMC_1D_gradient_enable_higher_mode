#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-project overrides for the CHT project (loaded automatically by
SetupData/main.py when it is run as:  python main.py CHT).

Only the names below differ from SetupData/parameters.py.
"""

# MC.C calls the Fortran forward model WITHOUT the OpenMP critical section
# (#pragma omp critical(forward) is commented out), so the threads share the
# Fortran COMMON blocks. On macOS that has been measured to crash outright most
# of the time at 4 threads, and it can corrupt forward results silently rather
# than crashing. 1 thread is slower but trustworthy.
# Raise this on CentOS, or once that race is fixed.
MC_number_of_cores = 1

# Higher mode (1st overtone) FORCED OFF for CHT, even though both run stations
# have a {sta}.hph file. This is the fundamental-mode-only control run; the
# higher-mode results it is compared against are archived in
# CHT/MonteCarlo_hph/. Set back to 1 to re-enable.
use_higher_mode = 0
