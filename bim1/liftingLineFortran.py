""""
    Question 2 from Homework 1:
        - Lifting line theory implementation in Fortran
"""
#%% Imports
import numpy as np
from llt import llt  # type: ignore

#%% Input test
twist = np.array([0., 0., 0., 0., 0., 0.])
gama = np.array([1., 1., 1., 1., 1., 1.])
span = 8.0
chord = 1.0
cl0 = 0.0
cla = 6.283185307179586
alpha = 0.08726646259971647
vinf = 10
rho_air = 1.225

res_llt, cl, cd = llt.llt_main(
    len(twist), 
    twist, 
    gama, 
    span, 
    chord, 
    cl0, 
    cla, 
    alpha, 
    vinf, 
    rho_air
    )

# print("Results from Fortran:")
# print("Lift coefficient (cl):", cl)
# print("Drag coefficient (cd):", cd) 
# print('Residues of the solution:', res_llt)