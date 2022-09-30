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
from modules.constants import *

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

### Transient linear simulations solver
# transientsolver = solver_linear (R_input, variable_input, bcs_input, dF_input)

def solver_linear (R_input, variable_input, bcs_input, dF_input): # ,form_compiler_parameters_input):
    problem = NonlinearVariationalProblem (R_input, variable_input, bcs = bcs_input, J = dF_input)#, form_compiler_parameters = form_compiler_parameters_input)

    transientsolver = NonlinearVariationalSolver (problem)

    prm = transientsolver.parameters
    # info (prm, True)

    # Nonlinear solver
    prm ['nonlinear_solver'] = 'newton'
    # prm ['print_matrix'] = False #True
    # prm ['print_rhs'] = False #True
    # prm ['symmetric'] = False #True

    # Linear solver
    # prm ['newton_solver']['absolute_tolerance'] = 1e-1 #1E-8
    # prm ['newton_solver']['convergence_criterion'] = 'residual' #'residual' 'incremental'
    # prm ['newton_solver']['error_on_nonconvergence'] = True
    prm ['newton_solver']['linear_solver'] = 'umfpack' # 'bicgstab' 'cg' 'gmres' 'minres' 'petsc' 'richardson' 'superlu_dist' 'tfqmr' 'umfpack'
    # prm ['newton_solver']['maximum_iterations'] = 10000
    # prm ['newton_solver']['preconditioner'] = 'ilu' # 'ilu' 'icc' 'petsc_amg' 'sor'
    # prm ['newton_solver']['relative_tolerance'] = 1e-1
    # prm ['newton_solver']['relaxation_parameter'] = 1.0
    # prm ['newton_solver']['report'] = True

    # Krylov solver
    # prm ['newton_solver']['krylov_solver']['absolute_tolerance'] = 1e-3 #1E-9
    # #     prm ['newton_solver']['krylov_solver']['error_on_nonconvergence'] = True
    # prm ['newton_solver']['krylov_solver']['maximum_iterations'] = 10000 # 500000
    # prm ['newton_solver']['krylov_solver']["monitor_convergence"] = True
    # prm ['newton_solver']['krylov_solver']["nonzero_initial_guess"] = True #False
    # prm ['newton_solver']['krylov_solver']['relative_tolerance'] = 1e-3
    # prm ['newton_solver']['krylov_solver']['report'] = True
        
    # LU solver
    # prm ['newton_solver']['lu_solver']['report'] = True
    # prm ['newton_solver']['lu_solver']['symmetric'] = False
    # prm ['newton_solver']['lu_solver']['verbose'] = True

    return transientsolver



### Transient Nonlinear simulations solver
def solver_nonlinear (R_input, variable_input, bcs_input, dF_input, form_compiler_parameters_input):
    problem = NonlinearVariationalProblem (R_input, variable_input, bcs = bcs_input, J = dF_input)#, form_compiler_parameters = form_compiler_parameters_input)

    transientsolver = NonlinearVariationalSolver (problem)
    prm = transientsolver.parameters
    # info (prm, True)

    # Nonlinear solver
    prm ['nonlinear_solver'] = 'newton'
    # prm ['print_matrix'] = False #True
    # prm ['print_rhs'] = False #True
    # prm ['symmetric'] = False #True

    # Linear solver
    # prm ['newton_solver']['absolute_tolerance'] = 1e-1 #1E-8
    # prm ['newton_solver']['convergence_criterion'] = 'residual' #'residual' 'incremental'
    # prm ['newton_solver']['error_on_nonconvergence'] = True
    prm ['newton_solver']['linear_solver'] = 'umfpack' # 'bicgstab' 'cg' 'gmres' 'minres' 'petsc' 'richardson' 'superlu_dist' 'tfqmr' 'umfpack'
    # prm ['newton_solver']['maximum_iterations'] = 10000
    # prm ['newton_solver']['preconditioner'] = 'ilu' # 'ilu' 'icc' 'petsc_amg' 'sor'
    # prm ['newton_solver']['relative_tolerance'] = 1e-1
    # prm ['newton_solver']['relaxation_parameter'] = 1.0
    # prm ['newton_solver']['report'] = True

    # Krylov solver
    # prm ['newton_solver']['krylov_solver']['absolute_tolerance'] = 1e-3 #1E-9
    # #     prm ['newton_solver']['krylov_solver']['error_on_nonconvergence'] = True
    # prm ['newton_solver']['krylov_solver']['maximum_iterations'] = 10000 # 500000
    # prm ['newton_solver']['krylov_solver']["monitor_convergence"] = True
    # prm ['newton_solver']['krylov_solver']["nonzero_initial_guess"] = True #False
    # prm ['newton_solver']['krylov_solver']['relative_tolerance'] = 1e-3
    # prm ['newton_solver']['krylov_solver']['report'] = True
        
    # LU solver
    # prm ['newton_solver']['lu_solver']['report'] = True
    # prm ['newton_solver']['lu_solver']['symmetric'] = False
    # prm ['newton_solver']['lu_solver']['verbose'] = True

    return prm ['nonlinear_solver'], prm ['newton_solver']['linear_solver']  


# (base) root@MacBook twofluidmodel # conda activate fenicsproject
# (fenicsproject) root@MacBook twofluidmodel # ./modules/transient_solvers.py
