# IMPORTS
import numpy as np
from fem_module import fem_module as fem # type: ignore
from fem_module_d import fem_module_d as fem_d # type: ignore
from fem_module_b import fem_module_b as fem_b # type: ignore

# EXECUTION

# Define test point
# Define areas
n_bars = 300
A0 = np.ones(n_bars)*0.5
d0 = np.linspace(0.0, 0.2, n_bars)**2

# FINITE-DIFFERENCE TEST

# Define step size
h = 1e-4

# Define perturbations of input variables
Ad = np.ones_like(A0)*0.1
dd = np.ones_like(d0)*0.02

# Compute reference point
res, C, V = fem.main(A0, d0)
# print(res)


# Compute perturbed point
res1, C1, V1 = fem.main(A0+h*Ad, d0+h*dd)

# Compute the directional derivatives with finite differences
resd_fd = (res1 - res)/h
Cd_fd = (C1 - C)/h
Vd_fd = (V1 - V)/h

# Use the Forward AD code to compute the same directional derivatives
res_tap, resd, C_tap, Cd, V_tap, Vd = fem_d.main_d(A0, Ad, d0, dd)

# Compare the norm of the difference among
# directional derivatives for each variable
FDtest_res = 1 - resd_fd/resd
FDtest_C = 1 - Cd_fd/Cd
FDtest_V = 1 - Vd_fd/Vd

print('FD test')
print('res:', np.max(np.abs(FDtest_res)))
print('C:', FDtest_C)
print('V:', FDtest_V)
print('')

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
resb_copy = resb.copy()
Cb_copy = Cb
Vb_copy = Vb

# Initialize input seeds for the derivative accumulation
Ab = np.zeros_like(A0)
db = np.zeros_like(d0)

# Execute the RAD code (remember to uses copies of output seeds)
fem_b.main_b(A0,Ab,d0,db,res,resb_copy,C,Cb_copy,V,Vb_copy) # Compute Ab e db

# Dot-product test
dotprod_inputs = 0.0
dotprod_inputs = dotprod_inputs + np.sum(Ad*Ab)
dotprod_inputs = dotprod_inputs + np.sum(dd*db)

dotprod_outputs = 0.0
dotprod_outputs = dotprod_outputs + np.sum(resd*resb)
dotprod_outputs = dotprod_outputs + Cd*Cb
dotprod_outputs = dotprod_outputs + Vd*Vb

dotprod_test = 1 - dotprod_inputs/dotprod_outputs

# Print results
print('DP test')
print('dotprod_inputs:',dotprod_inputs)
print('dotprod_outputs:',dotprod_outputs)
print('dotprod_test:',dotprod_test)
print('')
