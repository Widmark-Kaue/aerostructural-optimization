#%% Libs
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from cycler import cycler
from scipy.optimize import minimize

from utils import set_aiaa_style
from lltlib import LiftingLineOpt

#%% Plot settings
set_aiaa_style(16)

dpi = 600
format = 'pdf'
saveflag = True

#%% Path settings
rootdir = Path('.')
if not rootdir.absolute().name == 'llt_optimization':
    rootdir = Path('.','llt_optimization')


datadir = rootdir.joinpath('data')
datadir.mkdir(exist_ok=True, parents=True)
print("Data will be saved in:", datadir)

imagdir = rootdir.joinpath('images')
imagdir.mkdir(exist_ok=True)
print("Images will be saved in:", imagdir)

#%% Options for optmizer
tol = 1e-8
options = {'maxiter': 20000}

#%% Discretization settings and create LiftingLineOpt objects
vortices_list = [6, 8, 16,32,64]
lltclass_list = [LiftingLineOpt(nvortices=n) for n in vortices_list]

#%% Optimization - ITEM 2.5 and 2.6 -  6, 8 and 16 vortices
for i in range(0,len(lltclass_list)):
    lltclass = lltclass_list[i]
    
    # Constraints
    con = {
        'type': 'eq',
        'fun': lltclass.confun,
        'jac': lltclass.confungrad}

    # Bounds
    twist_min = np.deg2rad(-45)
    twist_max = np.deg2rad(45)
    bounds = [[twist_min,twist_max]]*lltclass.nvortices

    # Initial guess
    twist0 = np.zeros(lltclass.nvortices)

    result = minimize(lltclass.objfun,twist0,
                      jac = lltclass.objfungrad,
                      constraints=[con],
                      bounds=bounds,
                      method='SLSQP',
                      tol = tol,
                      options=options
                      )
    # save optimal twist
    lltclass.twist_opt = result.x
    
    print(f'------------------- {lltclass.nvortices} vortices -------------------')
    print(result)

    print('Optimal values:')
    print(f'twist = {np.round(np.rad2deg(result.x), 3)} [deg]')
    print(f'CD = {result.fun}')
    print(f'CL = {lltclass.solve_llt(result.x)[0]}')

#%% Compare results