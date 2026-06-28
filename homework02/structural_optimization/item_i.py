#%% Libs
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize
from strlib import StructuralOpt

from pathlib import Path
from utils import set_aiaa_style

#%% Plot settings
set_aiaa_style(16)

dpi = 600
format = 'pdf'
saveflag = True
#%% Data settings
saveHistory = False

#%% Path settings
rootdir = Path('.')
if not rootdir.absolute().name == 'structural_optimization':
    rootdir = Path('.','structural_optimization')


datadir = rootdir.joinpath('data')
datadir.mkdir(exist_ok=True, parents=True)
print("Data will be saved in:", datadir)

imagdir = rootdir.joinpath('images')
imagdir.mkdir(exist_ok=True)
print("Images will be saved in:", imagdir)

#%% Constants
E = 73.1e9          # [Pa]
sigma_y = 324e6     # [Pa]
rho = 2780          # [kg/m^3]
R = 0.05            # [m]
F1 = F2 = 1000      # [N]
rhoKS = 10          # 
ta = tb = 5         # [mm]

# Options for optmizer
tol = 1e-8
options = {'maxiter': 20000}

#%% Create strcutural optimize object
st = StructuralOpt(R,E,sigma_y,rho,F1,F2,rhoKS=rhoKS, save_history=True)

#%% Optimization - ITEM G

# Restrições
con = {
    'type': 'ineq',
    'fun': st.confun,
    'jac': st.confungrad}


# Bounds
bounds = [[0.1,50]]*2

# Initial guess
x0 = np.array([ta,tb])

result = minimize(st.objfun,x0, 
                jac=st.objfungrad,
                constraints=[con],
                bounds=bounds,
                method='slsqp',
                tol = tol,
                options = options)

opt_x_val = result.x
opt_m_val = result.fun
x_hist = np.array(st.x_hist)
f_hist = np.array(st.f_hist).reshape(-1,1)

print(result)
print()
print('Optimal valus:')
print(f' [ta, tb]: {np.round(result.x,4)} [mm]')
print(f'm = {result.fun:.4f} [kg]')
#%% Plot
gen = np.arange(1,len(f_hist)+1)

fig = plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.plot(gen,x_hist[:,0],'r', label = r'$t_a$')
plt.plot(gen,x_hist[:,1],'b',label = r'$t_b$')

plt.xlabel('Gen')
plt.ylabel(r'$\vec{x}$ [mm]')
plt.title('(a)')


plt.grid()
plt.legend()

plt.subplot(1,2,2)
plt.plot(gen,f_hist,'k')

plt.xlabel('Gen')
plt.ylabel(r'$m$ [kg]')
plt.title('(b)')
plt.grid()

plt.tight_layout()
plt.savefig(imagdir / f'q1_i.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
# plt.show()
# %%
