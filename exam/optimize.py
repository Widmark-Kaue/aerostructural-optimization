#%% libs
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize
from pathlib import Path

from asalib import ASAOptimization
from utils import set_aiaa_style

#%% Plot settings
set_aiaa_style(16)

dpi = 600
format = 'pdf'
saveflag = True

#%% Path settings
rootdir = Path('.')
if not rootdir.absolute().name == 'exam':
    rootdir = Path('.','exam')


datadir = rootdir.joinpath('data')
datadir.mkdir(exist_ok=True, parents=True)
print("Data will be saved in:", datadir)

imagdir = rootdir.joinpath('images')
imagdir.mkdir(exist_ok=True)
print("Images will be saved in:", imagdir)

#%% Wing geometry
Sref = 16.0 #[m^2]
AR = 6.0
b_wing = np.sqrt(AR*Sref)
cref = Sref/b_wing
npanels = 20

#%% Create ASA optimize object
asa = ASAOptimization(
    npanels=npanels,
    saveHistory=True, 
    span = b_wing,
    chords=cref,
    r = 0.05*cref
)


#%% Options for optmizer
tol = 1e-8
options = {'maxiter': 20000}

#%% Minimize Fuel burn

# clean history
asa.clean_history()

# Constraints
con1 = {'type': 'eq',
        'fun': asa.confun_eq,
        'jac': asa.confungrad_eq}

con2 = {'type': 'ineq',
        'fun': asa.confun_ineq,
        'jac': asa.confungrad_ineq}

con =  [con1, con2]


# Objective functions
objfun = lambda desVars: asa.objfun(desVars, func='FB')
objfungrad = lambda desVars: asa.objfungrad(desVars, func = 'FB')

# Initial guess 
twist0 = np.zeros(npanels, dtype=float)
t0  = np.ones(npanels,  dtype=float)*0.005
desVars0 = np.hstack([twist0,t0])
desVars0[npanels:] *= 100

# Bounds
twist_min = np.deg2rad(-45)
twist_max = np.deg2rad(45)
t_min = 0.001
t_max = asa.r

bounds = [[twist_min,twist_max]]*npanels + [[t_min*100,t_max*100]]*npanels

# Optimization
result = minimize(objfun, desVars0, 
                  jac = objfungrad,
                  constraints=con,
                  bounds=bounds,
                  method='SLSQP',
                  tol = tol,
                  options=options)

print(result)
# print(f"Optimal Fuel Burn: {result.fun * scale_obj:.4f} N")
#%% Plot history
twist_hist = asa.x_hist[:, :npanels]
t_hist = asa.x_hist[:,npanels:]
f_hist = asa.f_hist

plt.figure(figsize=(13,4))
plt.subplot(1,3,1)
plt.plot(twist_hist)

plt.xlabel(r'$N_{\text{fev}}$')
plt.ylabel(r'$\vec{\tau}$')

plt.grid()

plt.subplot(1,3,2)
plt.plot(t_hist*100)


plt.xlabel(r'$N_{\text{fev}}$')
plt.ylabel( r'$\vec{t} \times 10^2$')

plt.grid()

plt.subplot(1,3,3)
plt.semilogy(f_hist)

plt.xlabel(r'$N_{\text{fev}}$')
plt.ylabel(r'$FB$')

plt.grid()

plt.tight_layout()
plt.show()





