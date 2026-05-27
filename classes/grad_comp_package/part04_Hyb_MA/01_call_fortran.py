'''
01_call_fortran.py
This script interfaces with the Fortran code
'''

# IMPORTS
import numpy as np
from fem_module import fem_module as fem # type: ignore

# EXECUTION

# Define number of bars
n_bars = 3

# Define areas
A0 = np.ones(n_bars)*0.5

# Define displacements
d0 = np.ones(n_bars)*0.3

# Use the main interface
res, C, V = fem.main(A0, d0)
print('res',res)
print('C',C)
print('V',V)
print('')
