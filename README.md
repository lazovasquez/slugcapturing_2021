# Slug capturing with finite-elements (under construction)
- Contents
- Installation
- How to run the code

# Contents
The code includes 17 Jupyter notebooks:

## 0. Libraries

## 1. Constants

## 2. Functions 

## 3. Equations
The two fluid model equations

## 4. Test cases 
## 5. Stratified smooth

# Installation
The code requires [FEniCS project](https://fenicsproject.org) 2019.1.0 to run.

# How to run the code?

## 1. Stratified smooth flow calculation

## 2. Well-posedness and stability analysis of stratified smooth flow
Differential flow pattern maps, dispersion figures and eigenspectra

## 3. Stiffness analysis
Eigenspectra of spatial discretization schemes

## 4. Discrete flow pattern maps 

## 5. Validated Flow pattern maps 

## 6. Nonlinear simulations

## 5. Vefification of the code
Numerical solution of KY equations.





























#  Slug capturing model for gas-liquid flows in hilly terrain pipes
Model written in Python3 and C++ through the finite element packages [FEniCS](https://fenicsproject.org) and [Firedrake](https://www.firedrakeproject.org/).

## System Requirements
1. Python 3 → [Anaconda Distribution](https://www.anaconda.com/distribution/)
2. FEniCS 2019.1 → [FEniCS on Anaconda](https://fenicsproject.org/download/)
3. Firedrake → [Obtaining Firedrake](https://www.firedrakeproject.org/download.html)

## Codes 
- ### Linear stability analysis ###
1. Well-posedness of the two-fluid model (IKH)
2. Global stability analysis of the differential two-fluid model (IKH and VKH)
3. Fouier analysis of the differential two-fluid model (VKH)

- ### Taylor-Hood elements with high-order basis functions ###
4. Stiffness of the semi-discretized two-fluid model
5. Linear simulations of the fully discretized two-fluid model
6. Nonlinear simulations of the fully discretized two-fluid model

- ### Verification of the numerical method (Taylor-Hood elements) ###
7. Viscous Burguer's equation
8. Water faucet problem
9. Kreiss-Ystrom (KY) equations 
10. Incompressible two-fluid model

- ### Verification of the numerical method (Discontinuous Galerkin elements) ###
11. Viscous Burguer's equation
12. Water faucet problem
13. Kreiss-Ystrom (KY) equations 
14. Incompressible two-fluid model

- ### Verification of the numerical method (Spectral elements)###
15. Kreiss-Ystrom (KY) equations 

