'''
06_optimize.py
This script solves the optimization problem
'''

# IMPORTS
import numpy as np
from fem_module import fem_module as fem # type: ignore
from fem_module_b import fem_module_b as fem_b # type: ignore
from scipy.optimize import root, minimize
import time

# EXECUTION

# Define number of bars
n_bars = 40

# Define areas
A0 = np.ones(n_bars)*0.6

# Define displacements
d0 = np.ones(n_bars)*0.3

### OBJECTIVE FUNCTION

def objfun(A):

    # Define function that computes residuals for given displacements
    def resfunc(d):
        res, C, V = fem.main(A,d)
        return res

    # Solve the FEM problem
    sol = root(resfunc, d0)

    # Get solved displacements
    d = sol.x

    # Compute compliance function for the correct displacements
    res, C, V = fem.main(A,d)

    return C

def objfungrad(A):

    # Define function that computes residuals for given displacements
    def resfunc(d):
        res, C, V = fem.main(A,d)
        return res

    # Solve the FEM problem
    sol = root(resfunc, d0)

    # Get solved displacements
    d = sol.x

    # Initialize output variables for reverse calls
    # (the value does not matter since Fortran will overwrite them)
    C = 0.0
    V = 0.0
    res = np.zeros(n_bars)

    # Define function to compute the residual of the adjoint equation

    def adjfunc(resb):

        # Define seeds
        Ab = np.zeros_like(A)
        db = np.zeros_like(d)
        Cb = 1.0
        Vb = 0.0

        # Run reverse mode
        fem_b.main_b(A, Ab, d, db, res, resb.copy(), C, Cb, V, Vb)

        # Compute the residual of the adjoint equation
        adj_res = db

        return adj_res

    # Solve the adjoint problem
    psi0 = np.ones(n_bars)*1.0
    sol = root(adjfunc, psi0)

    # Get adjoint variables
    psi = sol.x

    # Initialize derivative seeds
    Ab = np.zeros_like(A)
    db = np.zeros_like(d)
    Cb = 1.0
    Vb = 0.0
    resb = psi

    # Run the reverse mode
    fem_b.main_b(A, Ab, d, db, res, resb, C, Cb, V, Vb)

    # Get total derivative
    dCdA = Ab

    return dCdA

### CONSTRAINT FUNCTION

def confun(A):

    # Define dummy values for displacement (since V does not depend on it)
    d = np.zeros_like(A)

    # Compute V
    res, C, V = fem.main(A,d)

    con = 1.0 - V

    return con

def confungrad(A):

    d = np.zeros_like(A)

    Ab = np.zeros(n_bars)
    db = np.zeros(n_bars)
    res = np.zeros_like(d)
    resb = np.zeros_like(d)
    C = 0.0 # Does not matter
    Cb = 0.0
    V = 0.0 # Does not matter
    Vb = 1.0
    fem_b.main_b(A, Ab, d, db, res, resb, C, Cb, V, Vb)
    dVdA = Ab

    congrad = -dVdA

    return congrad

### OPTIMIZATION

# Create list of constraints

con1 = {'type': 'ineq',
        'fun': confun,
        'jac': confungrad}

# Set lower bounds for areas
bounds = [[0.0, None]]*n_bars

# Run optimizer
result = minimize(objfun, A0, jac=objfungrad,
                  constraints=[con1],
                  bounds=bounds,
                  method='slsqp')

# Compute function values at the optimum
A = result.x
C = objfun(A)
V = 1.0-confun(A)

# Print results
print(result)
print('C:',C)
print('V:',V)
import matplotlib.pyplot as plt
fig = plt.figure()
plt.plot(result.x,'-o')
plt.show()
