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

if __name__ == '__main__':

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Transient simulator for two-phase flows')
    parser.add_argument('--case',
                        action='store',
                        metavar='case',
                        type=int,
                        help='Set of parameters.')
    arg = parser.parse_args()

    # Print Centerlines Modification Module
    print('')
    print('TWOFLUIDMODEL::SAYS:')

    # Load case parameters from case_x.json
    with open(os.path.join('cases', f'case_{arg.case}.json'), mode='r') as file1:
        case_data = json.load(file1)  # encoding='utf-8'

    # Setup
    simulation = case_data['setup']['simulation']

    # Equations parameters
    system = case_data['setup']['equations']['system']
    viscous_terms = case_data['setup']['equations']['viscous_terms']
    # Boundary conditions
    dirichlet_type = case_data['setup']['equations']['dirichlet_type']
    # Initial conditions
    IBVP = case_data['setup']['equations']['IBVP']
    effect = case_data['setup']['equations']['effect']

    # Plot reference conditions
    show_data = case_data['setup']['visualization']['show_data']

    # Phasic properties
    rho_l = case_data['phasic_properties']['liquid']['properties']['density']  # kg m^-3
    mu_l = case_data['phasic_properties']['liquid']['properties']['dynamic_viscosity']  # Pa s

    mu_g = 1.8e-5  # Pa s, gas viscosity# Phase properties
    c_g = case_data['phasic_properties']['gas']['properties']['compressibility']
    var4_0 = case_data['phasic_properties']['interface']['outlet_pressure']  # Pa

    # Numerical method
    elementspace = case_data['numerical_method']['discretization']['space']['elementspace']
    p = case_data['numerical_method']['discretization']['space']['order']
    CFL = case_data['numerical_method']['discretization']['CFL']
    nx = case_data['numerical_method']['discretization']['space']['elements_number']

    # Transient simulations
    transient_eigenspectrum = case_data['stability']['transient_spectrum']

    # Pipe inclination
    inclination = case_data['geometry']['inclination']

    if simulation == 2:
        discretization = 1
    elif any([simulation == 3, simulation == 4, simulation == 5]):
        discretization = 2

# (base) root@MacBook twofluidmodel # conda activate fenicsproject
# (fenicsproject) root@MacBook twofluidmodel # ./transient_simulator.py --case 1
