#%% libs
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

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
readData = True

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
asalistWeight = [ASAOptimization(
                    npanels=npanels,
                    saveHistory=True, 
                    span = b_wing[i],
                    chords=cref[i],
                    r = 0.05*cref[i]
                    )
     for i in range(len(AR))]

#%% Create list to store results
resultsWeight = []

#%% 5.7 Weight minimization

# Initial guess 
twist0 = np.zeros(npanels, dtype=float)
t0  = np.ones(npanels,  dtype=float)*0.005
desVars0 = np.hstack([twist0,t0]) # baseline
desVars0[npanels:] *= 100   

# Bounds
twist_min = np.deg2rad(-45)
twist_max = np.deg2rad(45)
t_min = 0.001
t_max = asalistWeight[0].r

bounds = [[twist_min,twist_max]]*npanels + [[t_min*100,t_max*100]]*npanels

# Loop for AR cases
for i in range(len(AR)):
    # get asa object
    asa = asalistWeight[i]
    
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
    objfun = lambda desVars: asa.objfun(desVars, func='Weight')
    objfungrad = lambda desVars: asa.objfungrad(desVars, func = 'Weight')
    
    
    # Optimization
    if not readData:
        result = minimize(objfun, desVars0, 
                          jac = objfungrad,
                          constraints=con,
                          bounds=bounds,
                          method='SLSQP')
        
        print(result)
        
        # Store results
        resultsWeight.append(deepcopy(result))
        
        # Export optimization data
        if exportData:
            res=deepcopy(result)
            res.asa=asa
            with open(datadir / f'Weight_results_AR{AR[i]}.pkl','wb') as f:
                dump(res,f)
            
            print(f'Weight minimization exported (AR={AR[i]})...')
    else:
        # Read optimization data
        with open(datadir / f'Weight_results_AR{AR[i]}.pkl', 'rb') as f:
            result = load(f)
        asalistWeight[i] = result.pop('asa')
        resultsWeight.append(result)
        print(f'Load Weight results(AR={AR[i]})...')
        
        
#%% Get opt values for AR = 6

# Get objects and results
asa = asalistWeight[0]
desVarsOpt = resultsWeight[0].x

twistOptdeg = np.rad2deg(desVarsOpt[:asa.npanels]) # convert to deg
tOptmm = desVarsOpt[asa.npanels:]/100 * 1000 # convert to mm

# Compute CL and Gamma from optimized wing
out, stateVars = asa.solve_asa(desVarsOpt)
CL = out['cl']
Gamma = stateVars[:npanels]

#%% 5.7 Baseline x Optimized
# close figures
plt.close('all')


twist0deg = np.rad2deg(twist0) # Convert to deg
t0mm = t0 * 1000 # Convert to mm

## Plot Design variables
plt.figure(figsize=(10,4))
# plot twist
plt.subplot(1,2,1)
plt.plot(asa.ypanels/b_wing[0], twistOptdeg,'bo-', label = 'Otimizado')
plt.plot(asa.ypanels/b_wing[0], twist0deg,'ko-', label = 'Baseline')
 
plt.xlim(-0.5,0.5)

plt.title('(a)')
plt.xlabel(r'$y/b$ [-]')
plt.ylabel(r'$\tau$ [deg]')


plt.grid()
plt.legend()

# plot t
plt.subplot(1,2,2)
plt.plot(asa.ypanels/b_wing[0], tOptmm,'bo-', label = 'Otimizado')
plt.plot(asa.ypanels/b_wing[0], t0mm,'ko-', label = 'Baseline')

plt.xlim(-0.5,0.5)

plt.title('(b)')
plt.xlabel(r'$y/b$ [-]')
plt.ylabel(r'$t$ [mm]')

plt.grid()
plt.legend()

plt.tight_layout()

plt.savefig(imagdir / f'q5.7_1_comp_desVars_baseline_x_opt.{format}', 
            dpi=dpi, bbox_inches='tight') if saveflag else None

# plt.show()


## Print Objfun and Constraints

obj0 = asa.objfun(desVars0, func = 'Weight')
objOPT = asa.objfun(desVarsOpt, func = 'Weight')

print('Objective function:')
print(f'Baseline: Weight = {obj0}')
print(f'Optimized: Weight = {objOPT}')

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

#%% 5.7 Comparison with eliptical lift distribution

# Compute eliptical lift distribution
span = asa.span
vinf = asa.Vinfm

y = asa.ypanels
gamma_0 = 2*vinf*Sref/(np.pi*span) * CL
gamma_elliptical = 2*gamma_0 * np.sqrt((0.5+ y/span) * (0.5-y/span))
# CD_an = cl**2*Sref/(np.pi*span**2)

plt.figure(figsize=(8,5))
plt.plot(y, gamma_elliptical, 'ko-', label = 'Elíptica')
plt.plot(asa.ypanels, Gamma, 'bo-', label = 'Otimizado')

plt.xlabel(r'$y/b$ [-]')
plt.ylabel(r'$\Gamma$ [m$^2$/s]')


plt.grid()
plt.legend()

plt.tight_layout()

plt.savefig(imagdir / f'q5.7_2_comp_opt_x_eliptical.{format}', 
            dpi=dpi, bbox_inches='tight') if saveflag else None
plt.show()

#%% 5.7 Plot margins of FEM nodes

margins = out['margins']
ymargins = asa.ymargins/asa.span

plt.figure(figsize=(8, 5))
plt.plot(ymargins, margins, 'ko')

for i in range(0,len(ymargins)+1,2):
    plt.plot(ymargins[i:i+2], margins[i:i+2])
    

plt.xlabel(r'$y/b$ [-]')
plt.ylabel(r'$m=1-\dfrac{\sigma}{\sigma_Y}$ [-]')


plt.grid()
plt.tight_layout()

plt.savefig(imagdir / f'q5.7_3_margins.{format}', 
            dpi=dpi, bbox_inches='tight') if saveflag else None
#%% 5.7 Weight Fractions

W0 = out['weight']
FB = out['fb']
Ws = W0 - FB - asa.fixedMass

weight_fraction = np.array([asa.fixedMass, Ws, FB])/W0 * 100
poslabel= [r'$W_f$',r'$W_s$', r'$W_{FB}$']
colors = ['b', 'r', 'g']

plt.figure(figsize=(8, 5))
for i in range(len(weight_fraction)):
    bars = plt.bar(poslabel[i],weight_fraction[i], color = colors[i],edgecolor =  'k')
    
    plt.bar_label(bars, fmt='%.1f%%', padding=3)
    
plt.bar([r'$W_0$'], 100, color='gray', edgecolor='k')

plt.xlabel('Weight')
plt.ylabel('Fração de Peso [%]')

# plt.grid()

plt.tight_layout()

plt.savefig(imagdir / f'q5.7_4_weight_fractions.{format}', 
            dpi=dpi, bbox_inches='tight') if saveflag else None
# plt.xticks(pos, labels=poslabel)

#%% plot history
plt.close('all')
# Get history
twist_hist = np.rad2deg(asa.x_hist[:,:npanels])
t_hist =  asa.x_hist[:, npanels:]*1000
weight_hist = asa.f_hist
liftExcess_hist = asa.h_hist
KSmargin_hist = asa.g_hist


plt.figure(figsize = (13,4))

N_total = len(weight_hist)
x1, x2 = max(0, N_total - 30), N_total +1

# Subplot 1: Twist
ax1 = plt.subplot(1,3,1)
plt.plot(twist_hist)

plt.title('(a)')
plt.xlabel(r'$N_{\text{fev}}$')
plt.ylabel(r'$\vec{\tau}$ [deg]')

plt.grid()

# Inset zoom for Twist
axins1 = inset_axes(ax1, width="40%", height="40%", loc="upper right")
axins1.plot(twist_hist)

y1 = np.min(twist_hist[x1:x2])
y2 = np.max(twist_hist[x1:x2])

margin_y = (y2 - y1) * 0.1 if y2 > y1 else 0.05
axins1.set_xlim(x1, x2)
axins1.set_ylim(y1 - margin_y, y2 + margin_y)

axins1.grid(True, linestyle=':', alpha=0.3)

axins1.tick_params(axis='both', which='major', labelsize=8)
mark_inset(ax1, axins1, loc1=2, loc2=4, fc="none", ec="0.5")

# Subplot 2: Thickness
ax2 = plt.subplot(1,3,2)
plt.plot(t_hist)

plt.title('(b)')
plt.xlabel(r'$N_{\text{fev}}$')
plt.ylabel(r'$\vec{t}$ [mm]')

plt.grid()

# Inset zoom for Thickness
axins2 = inset_axes(ax2, width="40%", height="40%", loc="upper right")
axins2.plot(t_hist)

y1 = np.min(t_hist[x1:x2])
y2 = np.max(t_hist[x1:x2])

margin_y = (y2 - y1) * 0.1 if y2 > y1 else 0.002
axins2.set_xlim(x1, x2)
axins2.set_ylim(y1 - margin_y, y2 + margin_y)

axins2.grid(True, linestyle=':', alpha=0.3)

axins2.tick_params(axis='both', which='major', labelsize=8)
mark_inset(ax2, axins2, loc1=2, loc2=4, fc="none", ec="0.5")

# Subplot 3: Weight
ax3 = plt.subplot(1,3,3)
plt.semilogy(weight_hist,  'k')

plt.title('(c)')
plt.xlabel(r'$N_{\text{fev}}$')
plt.ylabel(r'$W_0$ [N]')

plt.grid()

# Inset zoom for Weight (linear plot)
axins3 = inset_axes(ax3, width="40%", height="40%", loc="upper right")
axins3.plot(weight_hist, 'k')

y1 = np.min(weight_hist[x1:x2])
y2 = np.max(weight_hist[x1:x2])

margin_y = (y2 - y1) * 0.1 if y2 > y1 else 10.0
axins3.set_xlim(x1, x2)
axins3.set_ylim(y1 - margin_y, y2 + margin_y)

axins3.grid(True, linestyle=':', alpha=0.3)

axins3.tick_params(axis='both', which='major', labelsize=8)
mark_inset(ax3, axins3, loc1=2, loc2=4, fc="none", ec="0.5")

plt.tight_layout()

plt.savefig(imagdir / f'q5.7_5_desVars_obj_hist.{format}', 
            dpi=dpi, bbox_inches='tight') if saveflag else None
plt.show()


plt.figure(figsize=(9,4))
plt.subplot(1,2,1)
plt.semilogy(liftExcess_hist,'darkred')

plt.title('(a)')
plt.xlabel(r'$N_{\text{fev}}$')
plt.ylabel(r'$\Delta L$ [-]')

plt.grid()

plt.subplot(1,2,2)
plt.plot(KSmargin_hist,'darkgreen')

plt.title('(b)')
plt.xlabel(r'$N_{\text{fev}}$')
plt.ylabel(r'$m^{*}_{KS}$ [-]')

plt.grid()

plt.tight_layout()
plt.savefig(imagdir / f'q5.7_5_constraint_hist.{format}', 
            dpi=dpi, bbox_inches='tight') if saveflag else None

plt.show()
