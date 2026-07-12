#%% libs
import numpy as np
import matplotlib.pyplot as plt

from copy import deepcopy
from pickle import dump, load
from pathlib import Path
from scipy.optimize import minimize

from asalib import ASAOptimization
from utils import set_aiaa_style

#%% Plot settings
set_aiaa_style(16)

dpi = 600
format = 'pdf'
saveflag = True

#%% Data settings
exportData = True
readData = False

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
AR = np.array([6.0, 10])
b_wing = np.sqrt(AR*Sref)
cref = Sref/b_wing
npanels = 20

#%% Create ASA optimize object
asalistFB = [ASAOptimization(
                    npanels=npanels,
                    saveHistory=True, 
                    span = b_wing[i],
                    chords=cref[i],
                    r = 0.05*cref[i]
                    )
     for i in range(len(AR))]

#%% Create list to store results
resultsFB = []

#%% 5.6 Fuel burn minimization

# Initial guess 
twist0 = np.zeros(npanels, dtype=float)
t0  = np.ones(npanels,  dtype=float)*0.005
desVars0 = np.hstack([twist0,t0]) # baseline
desVars0[npanels:] *= 100   

# Bounds
twist_min = np.deg2rad(-45)
twist_max = np.deg2rad(45)
t_min = 0.001
t_max = asalistFB[0].r

bounds = [[twist_min,twist_max]]*npanels + [[t_min*100,t_max*100]]*npanels

# Loop for AR cases
for i in range(len(AR)):
    # get asa object
    asa = asalistFB[i]
    
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
    
    
    # Optimization
    if not readData:
        result = minimize(objfun, desVars0, 
                          jac = objfungrad,
                          constraints=con,
                          bounds=bounds,
                          method='SLSQP')
        
        print(result)
        
        # Store results
        resultsFB.append(deepcopy(result))
        
        # Export optimization data
        if exportData:
            res=deepcopy(result)
            res.asa=asa
            with open(datadir / f'FB_results_AR{AR[i]}.pkl','wb') as f:
                dump(res,f)
            
            print(f'FB minimization exported (AR={AR[i]})...')
    else:
        # Read optimization data
        with open(datadir / f'FB_results_AR{AR[i]}.pkl', 'rb') as f:
            result = load(f)
        asalistFB[i] = result.pop('asa')
        resultsFB.append(result)
        print(f'Load FB results(AR={AR[i]})...')

#%% 5.6 Baseline x Optimized
# close figures
plt.close()

# Get objects and results
asa = asalistFB[0]
desVarsOpt = resultsFB[0].x

twistOptdeg = np.rad2deg(desVarsOpt[:asa.npanels]) # convert to deg
tOptmm = desVarsOpt[asa.npanels:]/100 * 1000 # convert to mm

twist0deg = np.rad2deg(twist0) # Convert to deg
t0mm = t0 * 1000 # Convert to mm

## Plot Design variables
plt.figure(figsize=(10,4))
# plot twist
plt.subplot(1,2,1)
plt.plot(asa.ypanels/b_wing[0], twistOptdeg,'bo-', label = 'Optimized')
plt.plot(asa.ypanels/b_wing[0], twist0deg,'ko-', label = 'Baseline')

plt.xlim(-0.5,0.5)

plt.xlabel(r'$y/b$ [-]')
plt.ylabel(r'$\tau$ [deg]')


plt.grid()
plt.legend()

# plot t
plt.subplot(1,2,2)
plt.plot(asa.ypanels/b_wing[0], tOptmm,'bo-', label = 'Optimized')
plt.plot(asa.ypanels/b_wing[0], t0mm,'ko-', label = 'Baseline')

plt.xlim(-0.5,0.5)

plt.xlabel(r'$y/b$ [-]')
plt.ylabel(r'$t$ [mm]')

plt.grid()
plt.legend()

plt.tight_layout()

plt.savefig(imagdir / f'q5.6_1_comp_desVars_baseline_x_opt.{format}', 
            dpi=dpi, bbox_inches='tight') if saveflag else None

# plt.show()


## Print Objfun and Constraints

obj0 = asa.objfun(desVars0, func = 'FB')
objOPT = asa.objfun(desVarsOpt, func = 'FB')

print('Objective function:')
print(f'Baseline: FB = {obj0}')
print(f'Optimized: FB = {objOPT}')

deltaL0 = asa.confun_eq(desVars0)
deltaLOpt = asa.confun_eq(desVarsOpt)

KSmargin0 = asa.confun_ineq(desVars0)
KSmarginOpt = asa.confun_ineq(desVarsOpt)

print()
print('Constraints:')
print(f'Baseline: litfExcess = {deltaL0}')
print(f'Optimized: liftExcess = {deltaLOpt}')
print()
print(f'Baseline: KSmargin = {KSmargin0}')
print(f'Optimized: KSmargin = {KSmarginOpt}')


#%% plot history
twist_hist = asa.x_hist[:,:npanels]
t_hist =  asa.x_hist[:, npanels:]
fb_hist = asa.f_hist

plt.figure(figsize = (13,4))
plt.subplot(1,3,1)
plt.plot(twist_hist)
plt.ylabel(r'$\vec{\tau}$')

plt.subplot(1,3,2)
plt.plot(t_hist)
plt.ylabel(r'$\vec{t}$')

plt.subplot(1,3,3)
plt.semilogy(fb_hist)
plt.ylabel(r'$FB$')

plt.tight_layout()
plt.show()
