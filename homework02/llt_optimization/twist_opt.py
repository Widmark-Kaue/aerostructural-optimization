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
#%% Data settings
saveHistory = True

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


#%% Create LiftingLineOpt object
lltclass6 = LiftingLineOpt(nvortices=6, save_history=saveHistory)
lltclass8 = LiftingLineOpt(nvortices=8, save_history=saveHistory)
lltclass16 = LiftingLineOpt(nvortices=16, save_history=saveHistory)

cases = [lltclass6, lltclass8, lltclass16]

#%% Options for optmizer
# Options for optmizer
tol = 1e-8
options = {'maxiter': 20000}

#%% Optimization - ITEM 2.5 and 2.6 -  6, 8 and 16 vortices
for i in range(0,len(cases)):
    lltclass = cases[i]
    lltclass.clean_history()
    
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

    print(f'------------------- {lltclass.nvortices} vortices -------------------')
    print(result)

    print('Optimal values:')
    print(f'twist = {np.round(np.rad2deg(result.x), 3)} [deg]')
    print(f'CD = {result.fun}')
    print(f'CL = {lltclass.solve_llt(result.x)[0]}')

#%% Plot history
twist_hist = lltclass6.x_hist
CD_hist = lltclass6.f_hist

nfev = np.arange(1,len(CD_hist)+1)

fig = plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.plot(nfev,np.rad2deg(twist_hist),'-o')

plt.xlabel('nfev')
plt.ylabel(r'$\vec{\tau}$ [deg]')
plt.title('(a)')


plt.grid()
plt.legend([fr'$\tau_{i}$' for i in range(1,lltclass6.nvortices+1)],ncols = 2)

plt.subplot(1,2,2)
plt.plot(nfev,CD_hist,'ko-')

plt.xlabel('nfev')
plt.ylabel(r'$C_D$')
plt.title('(b)')
plt.grid()

plt.tight_layout()
plt.savefig(imagdir / f'q2_twist_opt.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
plt.show()