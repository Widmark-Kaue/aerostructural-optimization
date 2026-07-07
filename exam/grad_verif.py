#%% libs
import numpy as np
import matplotlib.pyplot as plt

from asa_module import asa_module as asa # type:ignore
from asa_module_d import asa_module_d as asa_d # type:ignore
from asa_module_b import asa_module_b as asa_b # type:ignore

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

#%% Seeds Foward method
# Input seeds
input_d = dict(
    twistd =  np.array([0.3, 0.1, 0.2, 0.4]),
    gamad = np.array([0.2, 0.3, 0.6, 0.4]),
    dd = np.ones_like(INPUTS['d'])*0.02,
    td = np.ones_like(INPUTS['t'])*0.1 
)
# Initizalize names for output seeds
outLabel_d = [f'{outLabel[i]}d' for i in range(7)]

# Call Foward AD method
outLabel_full = list(chain.from_iterable(zip(outLabel[:7], outLabel_d))) + outLabel[7:]
out_d = asa_d.asa_main_d(**INPUTS, **input_d)
out_d = dict(zip(outLabel_full,out_d))
out_d_slice = {kd: out_d[kd] for kd in outLabel_d}
pprint(out_d_slice)
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
    for kd in list(input_d.keys()):
        k = kd[:-1] # label without d
        INPUTS_cp[k] = INPUTS[k] + h[i] * input_d[kd]
    
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

#%% plot
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

plt.show()


# %%
