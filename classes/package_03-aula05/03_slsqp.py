'''
02_slsqp.py
Constrained Optimization
Cap. Ney Sêcco 2020
'''

# IMPORTS
import numpy as np
from scipy.optimize import minimize
import auxmod as am
import matplotlib.pyplot as plt

# EXECUTION
xlist = []
flist = []
g1list = []
g2list = []

# Define objective function

def objfun(x):
    f = 1-np.exp(-0.1*(x[0]**2 + x[1]**2))
    flist.append(f)
    return f

def objfungrad(x):
    gradf = np.array([np.exp(-0.1*(x[0]**2 + x[1]**2))*0.2*x[0],
                      np.exp(-0.1*(x[0]**2 + x[1]**2))*0.2*x[1]])
    return gradf

# Define constraints

def ineqcon(x):
    g1 = x[0]+x[1]-2
    g2 = x[1] - 1
    
    g1list.append(g1)
    g2list.append(g2)
    return [g1, g2]

def ineqcongrad(x):
    gradg1 = np.array([1,
                      1])
    gradg2 = np.array([0,
                       1])
    return [gradg1, gradg2]

def eqcon(x):
    h = x[0] - x[1]
    return h

def eqcongrad(x):
    gradh = np.array([1,
                      -1])
    return gradh

# Create list of constraints 
con1 = { 'type': 'ineq', 
        'fun': ineqcon,
        'jac': ineqcongrad}

con2 = { 'type': 'eq', 
        'fun': eqcon,
        'jac': eqcongrad}

cons = [con1, con2]

# Define starting point
x0 = np.array([3.0, 2.0])

# Define optimization bounds
bounds = [(-4,4),(-4,4)]

# Run optimizer
result = minimize(objfun, x0, jac = objfungrad, bounds=bounds, constraints=cons, method='SLSQP')

# Print results

print(result)                       
xopt = result.x

"""
    Nota-se que a segunda restrição é inativa visto na opção multipliers o segundo parâmetro é 0.
    Se retirarmos a segunda restrição da lista, o mínimo continua igual. 
"""

