'''
02_call_fortran.py
This script interfaces with the Fortran code
'''

# IMPORTS
import numpy as np
from fem_module import fem_module as fem # type: ignore
from fem_module_b import fem_module_b as fem_b # type: ignore
from scipy.optimize import root

# EXECUTION

# Define number of bars
n_bars = 3

# Define areas
A0 = np.ones(n_bars)*0.5

# Define displacements
d0 = np.ones(n_bars)*0.3


def resfunc(d):
    res, C, V = fem.main(A0, d)
    return res

sol = root(resfunc, d0)
print(sol)

d = sol.x

# Use the main interface
res, C, V = fem.main(A0, d)
print('res',res)
print('C',C)
print('V',V)
print('')

# ADJUTO C
def adjfunc_C(resb):
    Ab = np.zeros_like(A0)
    db = np.zeros_like(d)
    
    Cb = 1 
    Vb = 0
    
    fem_b.main_b(A0, Ab, d, db, res,resb.copy(), C, Cb, V, Vb)
    res_adj = db
    return db

psi0 = np.zeros_like(d) # variáveis adjuntas tem mesmo tamanho das variáveis de estado

sol = root(adjfunc_C, psi0)
print('ADJUNTO C')
print(sol)

psi = sol.x

Ab = np.zeros_like(A0)
db = np.zeros_like(d)    
Cb = 1 
Vb = 0

fem_b.main_b(A0, Ab, d, db, res, psi.copy(), C, Cb, V, Vb)

dCdA = Ab
print()
print('dCdA=', dCdA)

# ADJUTO V
def adjfunc_V(resb):
    Ab = np.zeros_like(A0)
    db = np.zeros_like(d)
    
    Cb = 0.0
    Vb = 1.0
    
    fem_b.main_b(A0, Ab, d, db, res,resb.copy(), C, Cb, V, Vb)
    res_adj = db
    return db

psi0 = np.zeros_like(d) # variáveis adjuntas tem mesmo tamanho das variáveis de estado

sol = root(adjfunc_V, psi0)
print('ADJUNTO V')
print(sol)

psi = sol.x

Ab = np.zeros_like(A0)
db = np.zeros_like(d)    
Cb = 0.0
Vb = 1.0

fem_b.main_b(A0, Ab, d, db, res, psi.copy(), C, Cb, V, Vb)

dVdA = Ab
print()
print('dVdA=', dVdA)