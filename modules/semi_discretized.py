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


## Semi-discretized
if simulation == 2:
    # Compute eigenspectra
    Real, Imag, A_stiffness, A_array = stiffness_function (Bm, visc, Cm, variable, dvariable)
    
    #################### Plot eigenvalues    
    plt.figure (5) # , figsize = mapsize
    fig, ax = plt.subplots ()
    ax.scatter (Real, 
                Imag, 
                s = area_scatter, 
                marker = listmarkers [0], 
                color = listcolor [4], 
                edgecolors = listcolor [0], 
                linewidths = line_width, 
                alpha = alphascatter)

    # ax.set_yscale ('symlog')

    # leg1 = ax.legend (loc = 'upper right', frameon = True, fontsize = 14);
    plt.grid (True, which = "both")
    plt.xlabel ('Re $(\mu)$ [s]', fontsize = label_size) # , fontsize = 18
    plt.ylabel ('Im $(\mu)$ [s]', fontsize = label_size) # , fontsize = 16
    # plt.xlim (-0.0025, -0.0015)
    # plt.ylim (-3, 3)
    matplotlib.rc ('xtick', labelsize = label_size)     
    matplotlib.rc ('ytick', labelsize = label_size)
    
    # Plot vertical line
    plt.axvline (0, label = 'pyplot vertical line', color = 'k')  

    # Save figure
    fig.set_size_inches (mapsize)
    plt.savefig('results/figures/semi_disc_maps/fig1.pdf',
                optimize = True,
                transparent = True,  
                dpi = dpi_elsevier)
    
    # Show figure
    plt.show
    
    # Information of eigenvalues
    print ("INFO: Max(Re[mu]) = ", max (Real))
    print ("INFO: Min(Re[mu]) = ", min (Real))
    print ("INFO: Max(Im[mu]) = ", max (Imag))
    print ("INFO: Min(Im[mu]) = ", min (Imag))

    #################### Plot acoustic waves
    plt.figure (6) # , figsize = mapsize
    fig, ax = plt.subplots ()
    ax.scatter (Real, 
                Imag, 
                s = area_scatter, 
                marker = listmarkers [0], 
                color = listcolor [4], 
                edgecolors = listcolor [0], 
                linewidths = line_width, 
                alpha = alphascatter)

    # plt.rcParams ['figure.figsize'] = mapsize
    # leg1 = ax.legend (loc = 'upper right', frameon = True, fontsize = 14);
    plt.grid (True, which = "both")
    plt.xlabel (r'Re ($\mu$) $[\it{s^{-1}}]$', fontsize = label_size)
    plt.ylabel (r'Im ($\mu$) $[\it{s^{-1}}]$', fontsize = label_size)
    plt.xlim (-0.08, -0.076)
    plt.ylim (-306.57, 306.57)
    matplotlib.rc ('xtick', labelsize = label_size)     
    matplotlib.rc ('ytick', labelsize = label_size)
    
    plt.axvline (0, label = 'pyplot vertical line', color = 'k')

    # Save figure
    fig.set_size_inches (mapsize)
    plt.savefig('results/figures/semi_disc_maps/fig2.pdf',
                optimize = True,
                transparent = True,  
                dpi = dpi_elsevier)
    
    # Show figure
    plt.show

    #################### Plot convective waves
    plt.figure (7)
    fig, ax = plt.subplots ()
    ax.scatter (Real, 
                Imag, 
                s = area_scatter, 
                marker = listmarkers [0], 
                color = listcolor [4], 
                edgecolors = listcolor [0], 
                linewidths = line_width, 
                alpha = alphascatter)

    # plt.rcParams ['figure.figsize'] = mapsize
    # leg1 = ax.legend (loc = 'upper right', frameon = True, fontsize = 14);
    plt.grid (True, which = "both")
    plt.xlabel (r'Re ($\mu$) $[\it{s^{-1}}]$', fontsize = label_size)
    plt.ylabel (r'Im ($\mu$) $[\it{s^{-1}}]$', fontsize = label_size)
    plt.xlim (-0.08, 0.02)
    plt.ylim (-30, 30)
    matplotlib.rc ('xtick', labelsize = label_size)     
    matplotlib.rc ('ytick', labelsize = label_size)
    
    plt.axvline (0, label = 'pyplot vertical line', color = 'k')

    # Save figure
    fig.set_size_inches (mapsize)
    plt.savefig('results/figures/semi_disc_maps/fig3.pdf',
                optimize = True,
                transparent = True,  
                dpi = dpi_elsevier)
    
    # Show figure
    plt.show

    #################### Stiffness matrix
    # Transform to numpy array
    A_mat = as_backend_type (A_stiffness).mat ()
    A_sparray = csr_matrix (A_mat.getValuesCSR ()[::-1], shape = A_mat.size)

    # Print matrix
    # print ("A_sparray = ", A_sparray)

    # Plot stiffness matrix
    plt.figure (8)
    fig, ax = plt.subplots ()
    plt.spy (A_sparray, color = 'k')
    plt.grid (True, which = "both")
    # # plt.xlim (0, l)
    # # ax.set_xlabel (r'L [m]', fontsize = 18)
    # # ax.set_ylabel (r'$p_i$ [Pa]', fontsize = 18)
    ax.xaxis.set_tick_params (which = 'major', direction = 'in', top = 'on')
    ax.xaxis.set_tick_params (which = 'minor', direction = 'in', top = 'on')
    ax.yaxis.set_tick_params (which = 'major', direction = 'in', right = 'on')
    ax.yaxis.set_tick_params (which = 'minor', direction = 'in', right = 'on')

    # Save figure
    fig.set_size_inches (mapsize)
    plt.savefig('results/figures/semi_disc_maps/fig4.pdf',
                optimize = True,
                transparent = True,  
                dpi = dpi_elsevier)
    
    # Show figure
    plt.show

# (base) root@MacBook twofluidmodel # conda activate fenicsproject
# (fenicsproject) root@MacBook twofluidmodel # ./modules/semi_discretized.py
