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

from scipy.stats import linregress
from pathlib import Path
from cycler import cycler

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
generateDataFlag = False                                    # if False, code will try read data from datadir 
#%%----------------------------------------------------------
### Path settings
rootdir = Path('.')
if not rootdir.absolute().name == 'homework01':
    rootdir = Path('.','homework01')

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
    

## Compute linear regression for log-log data
nm_linreg = linregress(np.log10(n), np.log10(nm_result[:, 1]))
de_linreg = linregress(np.log10(n[:-7]), np.log10(de_result[:-7, 1]))
cg_linreg = linregress(np.log10(n), np.log10(cg_result[:, 1]))
bfgs_linreg = linregress(np.log10(n), np.log10(bfgs_result[:, 1]))


## Build linear regression lines
n_fit = np.linspace(n[0], n[-1], 100)

nm_fit = np.polyval([nm_linreg.slope, nm_linreg.intercept], np.log10(n_fit)) # type: ignore
de_fit = np.polyval([de_linreg.slope, de_linreg.intercept], np.log10(n_fit)) # type: ignore
cg_fit = np.polyval([cg_linreg.slope, cg_linreg.intercept], np.log10(n_fit)) # type: ignore
bfgs_fit = np.polyval([bfgs_linreg.slope, bfgs_linreg.intercept], np.log10(n_fit)) # type: ignore

# Plot

# Dict with default parameters in context
colors = plt.color_sequences['Set1']
context = {'lines.linewidth': 2, 'lines.markersize': 6,'lines.linestyle': '--','lines.marker': 'o',
           'axes.prop_cycle': cycler('color', colors)}  # create cycle of color list

with plt.rc_context(context):
    plt.figure(figsize=(16, 4.5))
    
    plt.subplot(1,2, 1)
    
    plt.loglog(n, nm_result[:, 1], 'o',
               label = f'NM ($n=$${nm_linreg.slope:.2f}$, $R^2=$${nm_linreg.rvalue**2:.2f}$)') # type: ignore 
    plt.loglog(n, de_result[:, 1], 'o',
               label = f'DE ($n=$${de_linreg.slope:.2f}$, $R^2=$${de_linreg.rvalue**2:.2f}$)') # type: ignore
    plt.loglog(n, cg_result[:, 1], 'o',
               label = f'CG ($n=$${cg_linreg.slope:.2f}$, $R^2=$${cg_linreg.rvalue**2:.2f}$)') # type: ignore
    plt.loglog(n, bfgs_result[:, 1], 'o',
               label = f'BFGS ($n=$${bfgs_linreg.slope:.2f}$, $R^2=$${bfgs_linreg.rvalue**2:.2f}$)') # type: ignore

    # plot fit
    plt.plot(n_fit, 10**nm_fit, '-', color=colors[0])
    plt.plot(n_fit, 10**de_fit, '-', color=colors[1])
    plt.plot(n_fit, 10**cg_fit, '-', color=colors[2])
    plt.plot(n_fit, 10**bfgs_fit, '-', color=colors[3])
    
    plt.title('(a)')
    plt.xlabel('Number of dimensions')
    plt.ylabel('nfev')

    plt.legend()
    plt.grid(which='minor')
    
  
    
    # plt.figure(figsize=(8, 4.5))
    plt.subplot(1, 2, 2)
    plt.loglog(n, nm_result[:, 0], label = 'Nelder-Mead')
    plt.loglog(n, de_result[:, 0], label = 'Differential Evolution')
    plt.loglog(n, cg_result[:, 0], label = 'CG')
    plt.loglog(n, bfgs_result[:, 0], label = 'BFGS')
    
    plt.title('(b)')
    plt.xlabel('Number of dimensions')
    plt.ylabel(r'$|f(\mathbf{x}) - f^*|$')

    # plt.ylim([1e-20, 1])

    # plt.legend()
    
    plt.grid(which='both')
    
    plt.tight_layout()
    
    plt.savefig(imagdir / f'q3.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
    plt.show()

