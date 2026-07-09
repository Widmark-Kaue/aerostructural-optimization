#%% libs
import numpy as np
from asalib import ASAOptimization

from asa_module import asa_module as asa # type: ignore

#%% Print asa module
print(asa.asa_main.__doc__)

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
#%% Outputs
outname = ['resllt','resfem','liftexcess','margins','ksmargin','fb','weight','sref','cl']

# %%Call function
output = asa.asa_main(**geo, **aero,**struct,**oper)

outdic =dict(zip(outname,output))

keys = list(outdic.keys())

print('ASA MAIN TEST:')
for key in keys:
    print(f'{key} = {outdic[key]}\n')

# %% Validation of class
asaObj = ASAOptimization(4)
out = asaObj._run_asa(gama = aero['gama'], twist=aero['twist'], t=struct['t'], d = struct['d']) # type: ignore

print('ASA CLASS VALIDATION:')
print(30*'-')
for key in list(out.keys()):
    match = np.allclose(out[key], outdic[key])
    print(f'{key:<12}: Class = {out[key]}')
    print(f'{key:<12}: Call direct  = {outdic[key]}')
    print(f'Match {key:<10}: {match}\n')
