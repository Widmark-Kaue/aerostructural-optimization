#%% Libs
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from cycler import cycler
from scipy.stats import linregress
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

#%% funcs
def yc(lltclass:LiftingLineOpt):
    panel_length = lltclass.span/lltclass.nvortices
    idx = np.arange(1,lltclass.nvortices+1)
    yc  = 0.5*(panel_length - span) + (idx-1)*panel_length
    return yc

#%% Options for optmizer
tol = 1e-8
options = {'maxiter': 20000}

#%% Discretization settings and create LiftingLineOpt objects
vortices_list = np.array([6, 8, 16,32,64])
lltclass_list = [LiftingLineOpt(nvortices=n) for n in vortices_list]

# Variables for store optimal values
twist_opt_list = []
gama_opt_list = []
CL_opt_list = []
CD_opt_list = []

# Variables for store number of evaluations
nfev_list = []

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

    # Compute optimal CL and gama
    CL_opt, CD_opt, gama_opt = lltclass.solve_llt(result.x)
    
    # Store optimal values
    twist_opt_list.append(result.x)
    gama_opt_list.append(gama_opt)
    CL_opt_list.append(CL_opt)
    CD_opt_list.append(CD_opt)
    
    # Store number of evaluations
    nfev_list.append(result.nfev+result.njev) # sum of function and jacobian evaluations
    
    print(f'------------------- {lltclass.nvortices} vortices -------------------')
    print(result)

    print('Optimal values:')
    print(f'twist = {np.round(np.rad2deg(result.x), 3)} [deg]')
    print(f'CD = {result.fun}')
    print(f'CL = {CL_opt}')


# Convert to numpy arrays
CL_opt = np.array(CL_opt_list)
CD_opt = np.array(CD_opt_list)
nfev = np.array(nfev_list)

#%% Eliptical wing - Analytical solution
span = lltclass_list[0].span
chord = lltclass_list[0].chord
vinf = lltclass_list[0].Vinf
cl = lltclass_list[0].CL_target

Sref = span * chord

y = np.linspace(-span/2, span/2)
gamma_0 = 2*vinf*Sref/(np.pi*span) * cl
gamma_elliptical = 2*gamma_0 * np.sqrt((0.5+ y/span) * (0.5-y/span))
CD_an = cl**2*Sref/(np.pi*span**2)

print('Results from elliptical wing:')
print('Gamma values:', gamma_elliptical)
print("Lift coefficient (CL):", cl)
print("Drag coefficient (CD):", CD_an)


#%% Plot Circulation distribution
plt.figure(figsize=(8,5))
plt.plot(y,  gamma_elliptical,'k', label='Elliptical distribution')

for i in range(0,len(lltclass_list)):
    lltclass = lltclass_list[i]
    yc_i = yc(lltclass)
    plt.plot(yc_i, gama_opt_list[i], '-o',label=f'{lltclass.nvortices} vortices')


plt.xlabel('Spanwise location [m]')
plt.ylabel(r' $\Gamma$ [m$^2$/s]')

plt.legend(ncols =2)
plt.grid()
plt.tight_layout()
plt.savefig(imagdir / f'q2_comp_with_eliptical_circulation.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
# plt.show()

#%% Plot drag
Nv_inv = 1/vortices_list
opt_fit = linregress(Nv_inv, CD_opt)
Nv_fit = np.linspace(Nv_inv[0], 0, 100)

CD_opt_fit = np.polyval([opt_fit.slope, opt_fit.intercept], Nv_fit) #type: ignore
slope = opt_fit.slope #type: ignore
signal ='-' if np.sign(slope) == -1 else '+'

print('Fit results:')
print(f'CD = {opt_fit.intercept:.4f} {signal} {np.abs(opt_fit.slope):.4f} * Nv^{-1}') #type: ignore
print(f'R^2 = {opt_fit.rvalue**2:.4f}') #type: ignore


plt.figure(figsize=(14,5))
plt.subplot(1,2,1)
plt.plot(vortices_list, CD_opt, 'ko--', label='Optimization')
plt.axhline(CD_an, color='r', linestyle='--', label=f'Eliptical wing - $C_D = $ {CD_an:.4f}')

plt.xlabel('$N_v$')
plt.ylabel(r' $C_D$ ')
plt.title('(a)')

plt.legend()
plt.grid()

plt.subplot(1,2,2)
plt.plot(Nv_inv, CD_opt, 'ko', label='Optimization')
plt.plot(Nv_fit, CD_opt_fit, 'r-', label=f'Fit - $C_{{D,\\text{{opt}}}} = $ {opt_fit.intercept:.4f} {signal} {np.abs(opt_fit.slope):.4f} $N_v^{{-1}}$') #type: ignore

plt.xlabel('$N_v^{-1}$')
plt.ylabel(r' $C_{D,\text{opt}}$ ')
plt.title('(b)')

plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig(imagdir / f'q2_comp_with_eliptical_drag.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
# plt.show()

#%% Plot number of evaluations
plt.figure(figsize=(8,5))
plt.plot(vortices_list, nfev, 'ko--')

plt.xlabel('$N_v$')
plt.ylabel(r' $N_{\text{fev}}$ ')

plt.grid()

plt.tight_layout()
plt.savefig(imagdir / f'q2_comp_with_eliptical_nfev.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
plt.show()
