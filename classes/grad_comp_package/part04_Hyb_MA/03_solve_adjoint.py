'''
03_solve_adjoint.py
This script solves the FEM problem and solves the adjoint for C and V
'''

# IMPORTS
import numpy as np
from fem_module import fem_module as fem
from fem_module_b import fem_module_b as fem_b
from scipy.optimize import root

# EXECUTION

# Define number of bars
n_bars = 3

# Define areas
A0 = np.ones(n_bars)*0.5

# Define displacements
d0 = np.ones(n_bars)*0.3

# Define function that computes residuals for given displacements
def resfunc(d):
    res, C, V = fem.main(A0,d)
    return res

# Solve the FEM problem
sol = root(resfunc, d0)
print(sol)

# Get solved displacements
d = sol.x

# Compute function for the correct displacements
res,C,V = fem.main(A0,d)
print('res',res)
print('d',d)
print('C',C)
print('V',V)

### C ADJOINT

# Define function to compute the residual of the adjoint equation

def adjfunc(resb):

    # Define seeds
    Ab = np.zeros(n_bars)
    db = np.zeros(n_bars)
    Cb = ######
    Vb = ######

    # Run reverse mode
    fem_b.main_b(A0, Ab, d, db, res, resb.copy(), C, Cb, V, Vb)

    # Compute the residual of the adjoint equation
    adj_res = ######

    return adj_res

# Solve the adjoint problem
psi0 = np.ones(n_bars)*1.0
sol = root(adjfunc, psi0)
print('')
print(sol)

# Get adjoint variables
psi = sol.x

# Get total derivative
Ab = np.zeros(n_bars)
db = np.zeros(n_bars)
Cb = 1.0
Vb = 0.0
resb = psi.copy()
fem_b.main_b(A0, Ab, d, db, res, resb, C, Cb, V, Vb)

# Get total derivative
dCdA = Ab

print('')
print('dCdA',dCdA)

### V ADJOINT

# Define function to compute the residual of the adjoint equation

def adjfunc(resb):

    # Define seeds
    Ab = np.zeros(n_bars)
    db = np.zeros(n_bars)
    Cb = ######
    Vb = ######

    # Run reverse mode
    fem_b.main_b(A0, Ab, d, db, res, resb.copy(), C, Cb, V, Vb)

    # Compute the residual of the adjoint equation
    adj_res = ######

    return adj_res

# Solve the adjoint problem
psi0 = np.ones(n_bars)*1.0
sol = root(adjfunc, psi0)
print('')
print(sol)

# Get adjoint variables
psi = sol.x

# Get total derivative
Ab = np.zeros(n_bars)
db = np.zeros(n_bars)
Cb = 0.0
Vb = 1.0
resb = psi.copy()
fem_b.main_b(A0, Ab, d, db, res, resb, C, Cb, V, Vb)

# Get total derivative
dVdA = Ab

print('')
print('dVdA',dVdA)
