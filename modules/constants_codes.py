#!/usr/bin/env python3

from ast import Num
from fenics import *
from IPython.display import clear_output
from math import pi
from matplotlib import (rc, style)
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg as LA
from scipy import (linalg, matrix, sparse)
from scipy.interpolate import interp1d
from scipy.linalg import (eigvals, eig)
from scipy.misc import derivative as dtv
from scipy.optimize import (brenth, fsolve)
from scipy.sparse.linalg import eigs
from scipy.sparse import csr_matrix
from modules.constants_codes import *
from fenics import *
from scipy.misc import derivative as dtv

import math
import matplotlib
import matplotlib.font_manager
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pprint
import scipy.signal
import time
import ufl
import argparse
import json
import os
import subprocess
import h5py
import numpy as np


# Python parameters and functions
# Parameters for integration and functions for linearization.

# Gravity
g = 9.8  # Gravity

# Brenth interval guess
lima = DOLFIN_EPS
limb = 1 - lima

# Fsolve initial guess
x0 = DOLFIN_EPS  # 0.001


# For linearization of sources
def gradient(Cmat_element, nvariable):
    return dtv(Cmat_element, ref[nvariable - 1])


# 3D matrix for fourier analysis
def ThreeD(a, b, c):
    lst = [[[[] for col in range(a)] for col in range(b)] for row in range(c)]
    return lst


# For UFL
def Max(a, b):
    return (a + b + abs(a - b))/2


def Min(a, b):
    return (a + b - abs(a - b))/2

# FENICS PARAMETERS AND FUNCTIONS
# Parameters and functions

# Form compiler options
# parameters ['form_compiler']['representation'] = 'uflacs'
# parameters ["form_compiler"]["optimize"] = True
# parameters ["form_compiler"]["cpp_optimize"] = True
# parameters ["form_compiler"]['precision'] = 50
# parameters ["form_compiler"]["cpp_optimize_flags"] = "-O3 -ffast-math -march=native"
# parameters ["form_compiler"]["quadrature_degree"] = 8
# parameters ["allow_extrapolation"] = True
# parameters ["refinement_algorithm"] = "plaza_with_parent_facets"
# parameters ["std_out_all_processes"] = False; # Print log messages only from the root process in
# parallel

# For PETSc options 
# (https://fenicsproject.org/docs/dolfin/1.4.0/python/programmers-reference/cpp/la/
# SLEPcEigenSolver.html)
# PETScOptions.set ("st_ksp_type", "preonly")
# PETScOptions.set ("st_pc_type", "lu")
# PETScOptions.set ("st_pc_factor_mat_solver_package", "mumps")

# For linear and nonlinear solvers
# ffc_options = {"optimize": True, \
#                "eliminate_zeros": True, \
#                "precompute_basis_const": True, \
#                "precompute_ip_const": True}


# Test for PETSc and SLEPc
if not has_linear_algebra_backend("PETSc"):
    print("DOLFIN has not been configured with PETSc. Exiting.")
    exit()

if not has_slepc():
    print("DOLFIN has not been configured with SLEPc. Exiting.")
    exit()

# FILES
# Criation of folders and files.

# File for boundaries (linear simulations)
# file_boundaries = File ("domain_linear/boundaries.xml")
# # File for boundaries (nonlinear simulations)
# file_boundaries = File ("domain_nonlinear/boundaries.xml")# # File for nonlinear simulations


# File for linear simulations
ff_variable1 = File("results/fields/fields_linear/variable1.pvd", "compressed")
ff_variable2 = File("results/fields/fields_linear/variable2.pvd", "compressed")
ff_variable3 = File("results/fields/fields_linear/variable3.pvd", "compressed")
ff_variable4 = File("results/fields/fields_linear/variable4.pvd", "compressed")

# File for nonlinear simulations
ff_variable1_nonlinear = File("results/fields/fields_nonlinear/variable1.pvd", "compressed")
ff_variable2_nonlinear = File("results/fields/fields_nonlinear/variable2.pvd", "compressed")
ff_variable3_nonlinear = File("results/fields/fields_nonlinear/variable3.pvd", "compressed")
ff_variable4_nonlinear = File("results/fields/fields_nonlinear/variable4.pvd", "compressed")

# (base) root@MacBook twofluidmodel # conda activate fenicsproject
# (fenicsproject) root@MacBook twofluidmodel # ./modules/constants_codes.py
