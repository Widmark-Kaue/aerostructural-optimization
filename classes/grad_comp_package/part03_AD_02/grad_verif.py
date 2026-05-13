# IMPORTS
import numpy as np
from fem_module import fem_module as fem
from fem_module_d import fem_module_d as fem_d
from fem_module_b import fem_module_b as fem_b

# EXECUTION

# Define test point
# Define areas
n_bars = 3
A0 = np.ones(n_bars)*0.5
d0 = np.linspace(0.0, 0.2, n_bars)**2

# FINITE-DIFFERENCE TEST

# Define step size
h = ???

# Define perturbations of input variables
???

# Compute reference point
???

quit()

# Compute perturbed point
???

# Compute the directional derivatives with finite differences
resd_fd = ???
Cd_fd = ???
Vd_fd = ???

# Use the Forward AD code to compute the same directional derivatives
???

# Compare the norm of the difference among
# directional derivatives for each variable
FDtest_res = ???
FDtest_C = ???
FDtest_V = ???

print('FD test')
print('res:', FDtest_res)
print('C:', FDtest_C)
print('V:', FDtest_V)
print('')

quit()

# DOT-PRODUCT TEST

# Define random inputs seeds for FAD
Ad = np.ones_like(A0)*0.1
dd = np.ones_like(A0)*0.02

# Execute the FAD code
res, resd, C, Cd, V, Vd = fem_d.main_d(A0, Ad, d0, dd)

# Define random output seeds for RAD
resb = np.ones_like(res)*0.5
Cb = 0.3
Vb = 0.1

# Make copies of the output seeds because the RAD
# code will modify them, and this will break the
# dot-product test that we will do later on
???

# Initialize input seeds for the derivative accumulation
???

# Execute the RAD code (remember to uses copies of output seeds)
???

# Dot-product test
dotprod_inputs = 0.0
dotprod_inputs = dotprod_inputs + ???
dotprod_inputs = dotprod_inputs + ???

dotprod_outputs = 0.0
dotprod_outputs = dotprod_outputs + ???
dotprod_outputs = dotprod_outputs + ???
dotprod_outputs = dotprod_outputs + ???

dotprod_test = ???

# Print results
print('DP test')
print('dotprod_inputs:',dotprod_inputs)
print('dotprod_outputs:',dotprod_outputs)
print('dotprod_test:',dotprod_test)
print('')