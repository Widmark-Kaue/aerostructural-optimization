#%% libs
import numpy as np
import matplotlib.pyplot as plt

from copy import deepcopy

from utils import set_aiaa_style

from llt import llt # type: ignore - Código original
from llt_d import llt_d # type: ignore - Código direto
from llt_b import llt_b # type: ignore - Código reverso

#%% plot settings
set_aiaa_style(16)

#%% Inputs
inputs = dict(
            twist = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            gama = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
            span = 8.0,
            chord = 1.0,
            cl0 = 0.0,
            cla = 2*np.pi,
            alpha = 5*np.pi/180,
            vinf = 10,
            rho_air = 1.225,
            )


inputs_seed = dict(
                twistd = np.array([0.3, 0.1, 0.2, 0.4, 0.6, 0.5]),
                gamad = np.array([0.2, 0.3, 0.6, 0.4, 0.5, 0.1])
                )

#%% Finite Diference test
exp = np.arange(2,-12,-1,  dtype=float)
h = 10**(exp)
inputs1 = deepcopy(inputs)

FDtest_CL = np.zeros(len(exp))
FDtest_CD = np.zeros(len(exp))
FDtest_res = np.zeros((len(exp), len(inputs['twist']))) # type: ignore

# Compute gradient using AD Foward
_,res_lltd,_,CLd,_,CDd = llt_d.llt_main_d(**inputs,**inputs_seed)
    
print(f'CLd: {CLd}')
print(f'CDd: {CDd}')
print(f'res_lltd: {res_lltd}')
print()
for i in range(len(exp)): 
    inputs1['twist'] = inputs['twist'] + h[i]*inputs_seed['twistd']
    inputs1['gama'] = inputs['gama'] + h[i]*inputs_seed['gamad']

    # Compute gradient using FD
    res_llt,CL,CD= llt.llt_main(**inputs)
    res_llt1,CL1,CD1= llt.llt_main(**inputs1)

    res_lltd_FD = (res_llt1-res_llt)/h[i]
    CLd_FD = (CL1-CL)/h[i]
    CDd_FD = (CD1-CD)/h[i]
    
    # Finite Difference test
    FDtest_CL[i] = abs(1 - CLd_FD/CLd)
    FDtest_CD[i] = abs(1 - CDd_FD/CDd)
    FDtest_res[i] =np.abs(1 - res_lltd_FD/res_lltd
)    
    print(f'-----------h = {h[i]} ------------')
    print(f'CLd_FD: {CLd_FD}')
    print(f'CDd_FD: {CDd_FD}')
    print(f'res_lltd_FD: {res_lltd_FD}')
    print()
    print(f'FDtest_CL: {FDtest_CL[i]}')
    print(f'FDtest_CD: {FDtest_CD[i]}')
    print(f'FDtes_res: {FDtest_res[i]}')
    
    
plt.figure()
plt.subplot(1,2,1)
plt.loglog(h,FDtest_CL)
plt.loglog(h,FDtest_CD)
plt.xlim(max(h), min(h))

plt.grid()
plt.subplot(1,2,2)
plt.loglog(h,FDtest_res)
plt.xlim(max(h), min(h))

plt.grid()
plt.show()
# %%


# %%
