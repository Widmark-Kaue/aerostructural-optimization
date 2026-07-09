#%% libs
import numpy as np
import matplotlib.pyplot as plt

from asa_module import asa_module as asa # type:ignore
from asa_module_d import asa_module_d as asa_d # type:ignore
from asa_module_b import asa_module_b as asa_b # type:ignore
from asalib import ASAOptimization

from pprint import pprint
from itertools import chain
from copy import deepcopy
from utils import set_aiaa_style
from pathlib import Path

#%% plot settings
set_aiaa_style(16)
dpi = 600
format = 'pdf'
saveflag = True
#%% Path settings
rootdir = Path('.')
if not rootdir.absolute().name == 'exam':
    rootdir = Path('.','exam')

imagdir = rootdir.joinpath('images')
imagdir.mkdir(exist_ok=True)
print("Images will be saved in:", imagdir)

#%% Print funcs
print('ORIGINAL CODE')
print(100*'-')
print(asa.asa_main.__doc__)
print()
print('FOWARD AD')
print(100*'-')
print(asa_d.asa_main_d.__doc__)
print()
print('REVERSE AD')
print(100*'-')
print(asa_b.asa_main_b.__doc__)

#%% inputs
Vinfm = 60 # [m/s]
alpha = np.deg2rad(5) # [rad]

#Geometric parameters
geo = dict(
    span = 8 # [m]
)
# Aerodynamic parameters
aero =  dict(
    gama = np.array([80.0,80.0,80.0,80.0]),
    twist = np.array([0.0,0.0,0.0,0.0]),
    chords = np.array([1.0,1.0,1.0,1.0]),                           #[m]
    cl0 = np.array([0.0,0.0,0.0,0.0]), 
    cla = np.array([6.283,6.283,6.283,6.283]),                      #[1/rad]
    vinf= np.array([Vinfm*np.cos(alpha), 0, Vinfm*np.sin(alpha)]),  #[m/s]
    rhoair = 1.225                                                  #[kg/m^3]
)

# Structural parameters
struct = dict(
    r = np.array([0.1,0.1,0.1,0.1]),                                #[m]
    t = np.array([0.005,0.005,0.005,0.005]),                        #[m]
    e = 73.1e9,                                                       #[Pa]
    rhomat = 2780.0,                                                #[kg/m^3]
    sigmay = 324e6,                                                 #[Pa]
    pks = 200,
    d = np.array([0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001])  
)

# Operational parameters
oper = dict(
    cd0 = 0.0270,
    fixedmass = 700.0,                                              # [kg]
    g = 9.8,                                                        # [m/s2]
    endurance = 4.0*60.0*60.0,                                      # [s]
    tsfc = 0.5/3600.0,                                              # [1/s]
    loadfactor = 3.0*1.5
)

# Merge - Create input dict
INPUTS = geo | aero | struct | oper
#%% Outputs Label
outLabel = ['resllt','resfem','liftexcess','margins','ksmargin','fb','weight','sref','cl']
outLatex = ['r_{llt}', 'r_{fem}', r'\Delta L','m','m_{KS}','FB','W']
outLatex = dict(zip(outLabel,  outLatex))

# names for output seeds Foward AD
outLabel_d = [f'{outLabel[i]}d' for i in range(7)]

# names for output seeds Reverse AD
outLabel_b = [f'{outLabel[i]}b' for i in range(7)]

#%% Seeds Foward method
# Input seeds
inp_d = dict(
    twistd =  np.array([0.3, 0.1, 0.2, 0.4]),
    gamad = np.array([0.2, 0.3, 0.6, 0.4]),
    dd = np.ones_like(INPUTS['d'])*0.02,
    td = np.ones_like(INPUTS['t'])*0.1 
)

# Call Foward AD method
outLabel_full = list(chain.from_iterable(zip(outLabel[:7], outLabel_d))) + outLabel[7:]
out_d = asa_d.asa_main_d(**INPUTS, **inp_d)
out_d = dict(zip(outLabel_full,out_d))
out_d= {kd: out_d[kd] for kd in outLabel_d}
pprint(out_d)
# %% Finite Diferences Test

# Compute values for non-pertubation case
out0 = asa.asa_main(**INPUTS)
out0 = dict(zip(outLabel, out0)) 

# Temp inputs dict
INPUTS_cp = deepcopy(INPUTS)

# FD Steps
exp = np.arange(2,-21,-1,  dtype=float)
h = 10**(exp)

# Cread FDtest dict
FDtest = {kd:[] for kd in outLabel_d}

# For loop on h
for i in range(len(exp)):
    # Apply pertubation on INPUT variables
    for kd in list(inp_d.keys()):
        k = kd[:-1] # label without d
        INPUTS_cp[k] = INPUTS[k] + h[i] * inp_d[kd]
    
    # Compute pertubation values
    outFD = asa.asa_main(**INPUTS_cp)
    outFD = dict(zip(outLabel, outFD))
    
    # Compute FD test
    print(f'-----------h = {h[i]} ------------')
    for kd in outLabel_d:
        k = kd[:-1] # label without d
        FD = (outFD[k] - out0[k])/h[i]
        FDtest[kd].append(np.abs(1 - FD/out_d[kd]))
        
        print(f'{kd}: {FDtest[kd][-1]}')

#%% plot FDtest
plt.close()
plt.figure(figsize = (13,9))
hchoose = 10**(-8)
i = 1
for kd in outLabel_d:
    k = kd[:-1]
    fdtest = np.array(FDtest[kd])
    if fdtest.ndim == 1:
        plt.subplot(2,2,1)
        plt.loglog(h,fdtest, label = fr'${outLatex[k]}$')

        plt.axvline(hchoose,color = 'k', linestyle = '--')
        
        plt.xlabel('h [-]')
        plt.ylabel('FD test [-]')

        plt.xlim(max(h), min(h))
        plt.grid(True)
        plt.legend()
    else:
        i = i+1
        plt.subplot(2,2,i)
        plt.loglog(h, fdtest)
        plt.axvline(hchoose,color = 'k', linestyle = '--')
        
        plt.xlabel('h [-]')
        plt.ylabel('FD test [-]')

        plt.xlim(max(h), min(h))

        plt.grid()
        if outLatex[k].endswith('}'):
            label=[fr'${outLatex[k][:-1]},{i}{'}'}$' for i in range(1,len(out0[k])+1)]
        else:
            label = [fr'${outLatex[k]}_{i}$' for i in range(1,len(out0[k])+1)]
        ncols = 2 
       
        
        plt.legend(label,ncols= ncols, loc = 'best')
        
        
plt.tight_layout()
plt.savefig(imagdir / f'FDtest.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
# plt.show()

#%% Seeds for reverse method
rng = np.random.seed(42)
# Original outputs
out_reverse = {k: np.zeros_like(out0[k]) for k in outLabel}
print()
print('out_revers:')
pprint(out_reverse)

# Output seeds
out_b = []
for kb in outLabel_b:
    var = np.array(out_reverse[kb[:-1]])
    ndim = var.ndim
    
    if ndim == 0:
        out_b.append(np.random.rand(1)[0])
    else:
        out_b.append(np.random.rand(len(var)))

out_b = dict(zip(outLabel_b,out_b))
out_b_copy = deepcopy(out_b) # Make copies of the output seeds


print()
print('out_b:')
pprint(out_b)
    
# Initialize input seeds for derivative accumulation
inp_b = dict(
    twistb = np.zeros_like(INPUTS['twist']),
    gamab = np.zeros_like(INPUTS['gama']),
    db = np.zeros_like(INPUTS['d']),
    tb = np.zeros_like(INPUTS['t'])
)
print()
print('Before call reverse conde - inp_b:')
pprint(inp_b)

# Call reverse code
asa_b.asa_main_b(**INPUTS,**inp_b,**out_reverse,**out_b_copy)

print()
print('After call reverse conde - inp_b:')
pprint(inp_b)
# %% Dot product test
xd = np.concatenate(tuple(inp_d[kd] for kd in inp_d.keys()))
xb = np.concatenate(tuple(inp_b[kb] for kb in inp_b.keys()))

yd = np.concatenate(tuple(np.array(out_d[kd], ndmin = 1) for kd in out_d.keys()))
yb = np.concatenate(tuple(np.array(out_b[kb], ndmin = 1) for kb in out_b.keys()))

dotprod_inputs = np.sum(xd*xb)
dotprod_outputs = np.sum(yd*yb)

dotprod_test = 1 - dotprod_inputs/dotprod_outputs
print()
print('DP test')
print('dotprod_inputs:',dotprod_inputs)
print('dotprod_outputs:',dotprod_outputs)
print('dotprod_test:',dotprod_test)
print('')

# %% Call reverse method via ASAOptimization wrapper to validate it
opt = ASAOptimization(npanels=len(INPUTS['gama'])) # type: ignore

grad_wrapper = opt._run_asa_b(
    gama=INPUTS['gama'], # type: ignore
    twist=INPUTS['twist'], # type: ignore
    t=INPUTS['t'], # type: ignore
    d=INPUTS['d'], # type: ignore
    output_seeds=out_b
)

print('ASA REVERSE CLASS VALIDATION:')
print(30*'-')
for key in ['gamab', 'twistb', 'tb', 'db']:
    match = np.allclose(grad_wrapper[key], inp_b[key])
    print(f'{key:<8}: Class = {grad_wrapper[key]}')
    print(f'{key:<8}: Call direct  = {inp_b[key]}')
    print(f'Match {key:<6}: {match}\n')