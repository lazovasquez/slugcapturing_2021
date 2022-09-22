#!/usr/bin/env python3

from fenics import *
from scipy.misc import derivative as dtv


# Gravity
g = 9.8  # Gravity

transient_eigenspectrum = 1

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
