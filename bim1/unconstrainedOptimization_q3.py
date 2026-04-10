"""
    Author: Widmark Kauê Silva Cardoso
    
    Question 3 from homework 1: 
        - Unconstrained optimization of the Rosenbrock function
        
"""

#%%----------------------------------------------------------
### Imports
import os
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from hw01_q3_mod import nm_opt, de_opt, cg_opt, bfgs_opt
from utils import set_aiaa_style
#%%----------------------------------------------------------
### Clear screen
os.system('clear') # for linux

#%%----------------------------------------------------------
### plot settings
set_aiaa_style(16)

dpi = 600
format = 'pdf'
saveflag = True
#%%----------------------------------------------------------
### Data settings
generateDataFlag = True
#%%----------------------------------------------------------
### Path settings
rootdir = Path('.')
if not rootdir.absolute().name == 'bim1':
    rootdir = Path('.','bim1')

datadir = rootdir.joinpath('data')
datadir.mkdir(exist_ok=True, parents=True)
print("Data will be saved in:", datadir)

imagdir = rootdir.joinpath('images')
imagdir.mkdir(exist_ok=True)
print("Images will be saved in:", imagdir)
#%%----------------------------------------------------------
### 3.4 - Problem scale 
n = np.arange(2, 36)
tol = 1e-6

if generateDataFlag:
    nm_result = np.zeros((len(n), 2))
    de_result = np.zeros_like(nm_result)
    cg_result = np.zeros_like(nm_result)
    bfgs_result = np.zeros_like(nm_result)

    for i in range(len(n)):
        print(f'------ Number of design variables: {n[i]} ------')
        nm_result[i] = nm_opt(n[i], tol = tol)[1:] 
        de_result[i] = de_opt(n[i], tol = tol)[1:] 
        cg_result[i] = cg_opt(n[i], tol = tol)[1:] 
        bfgs_result[i] = bfgs_opt(n[i], tol = tol)[1:]


    np.savetxt(datadir / 'nm_result.dat', nm_result, header='f_opt\tnfev')
    np.savetxt(datadir / 'de_result.dat', de_result, header='f_opt\tnfev')
    np.savetxt(datadir / 'cg_result.dat', cg_result, header='f_opt\tnfev')
    np.savetxt(datadir / 'bfgs_result.dat', bfgs_result, header='f_opt\tnfev')
else:
    nm_result = np.loadtxt(datadir / 'nm_result.dat')
    de_result = np.loadtxt(datadir / 'de_result.dat')
    cg_result = np.loadtxt(datadir / 'cg_result.dat')
    bfgs_result = np.loadtxt(datadir / 'bfgs_result.dat')
    




# Function evaluation

