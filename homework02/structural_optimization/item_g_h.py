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
saveHistory = True

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


#%% Optimization - ITEM G
# Create strcutural optimize object
st = StructuralOpt(R,E,sigma_y,rho,F1,F2,rhoKS=rhoKS, save_history=True)

# Restrições
con = {
    'type': 'ineq',
    'fun': st.confunKS,
    'jac': st.confunKSgrad}

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

if saveHistory:
    data = np.concatenate((x_hist,f_hist), axis=1)
    np.savetxt(datadir / 'opt_KS_history.dat', data, header='ta\ttb\tm')
#%%
gis = st.confunKS(opt_x_val) 

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
plt.savefig(imagdir / f'q1_g.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
# plt.show()

#%% OPTIMIZATION - ITEM H
st.save_history = False

# rhoKSValues1 = np.arange(1.4,10,0.2)
# rhoKSValues2 = np.arange(10,101)
# rhoKSValues = np.concatenate((rhoKSValues1,rhoKSValues2))

rhoKSValues = np.logspace(0.15,2)
X = np.zeros((len(rhoKSValues),2))
M = np.zeros(len(rhoKSValues))

for i in range(len(rhoKSValues)):
    st.rhoKS = rhoKSValues[i]
    result = minimize(st.objfun,x0, 
                jac=st.objfungrad,
                constraints=[con],
                bounds=bounds,
                method='slsqp',
                tol = tol,
                options = options)
    
    print('---------------------------')
    print(f'rhoKS = {rhoKSValues[i]}')
    print(result.message)
    if result.success:
        X[i] = result.x
        M[i] = result.fun
#%% plot
X_mod = np.zeros_like(X)
X_mod[:,0] =  X[:,0]/opt_x_val[0]
X_mod[:,1] =  X[:,1]/opt_x_val[1]

M_mod =  M/opt_m_val

fig = plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
# plt.plot(1/rhoKSValues,X[:,0]/opt_x_val[0],'ro', label = r'$t_a/t_{a,\text{ref}}$')
# plt.plot(1/rhoKSValues,X[:,1]/opt_x_val[1],'bo',label = r'$t_b/t_{b,\text{ref}}$')
plt.semilogx(1/rhoKSValues,X_mod[:,0],'ro', label = r'$t_a/t_{a,\text{ref}}$')
plt.semilogx(1/rhoKSValues,X_mod[:,1],'bo',label = r'$t_b/t_{b,\text{ref}}$')

plt.xlabel(r'$1/\rho_{KS}$')
plt.ylabel(r'$\vec{x}/\vec{x}_{\text{ref}}$')
plt.title('(a)')

plt.ylim(0,9)
plt.yticks(np.arange(0,9))
plt.grid()
plt.legend()

plt.subplot(1,2,2)
# plt.plot(1/rhoKSValues,M/opt_m_val,'ko')
plt.semilogx(1/rhoKSValues,M_mod,'ko')

plt.xlabel(r'$1/\rho_{KS}$')
plt.ylabel(r'$m/m_{\text{ref}}$')
plt.title('(b)')

plt.ylim(0,7)

plt.grid()

plt.tight_layout()
plt.savefig(imagdir / f'q1_h.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
plt.show()
# %%
