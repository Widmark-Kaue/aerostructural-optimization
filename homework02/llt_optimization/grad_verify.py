#%% libs
import numpy as np
import matplotlib.pyplot as plt

from copy import deepcopy
from pathlib import Path

from utils import set_aiaa_style

from llt import llt # type: ignore - Código original
from llt_d import llt_d # type: ignore - Código direto
from llt_b import llt_b # type: ignore - Código reverso

#%% plot settings
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

#%% Call origional func
res_llt_OR, CL_OR, CD_OR = llt.llt_main(**inputs)
print('Original function:')
print(f'res: {res_llt_OR}')
print(f'CL: {CL_OR}')
print(f'CD: {CD_OR}')

#%% Finite Diference test
inputs_seed = dict(
                twistd = np.array([0.3, 0.1, 0.2, 0.4, 0.6, 0.5]),
                gamad = np.array([0.2, 0.3, 0.6, 0.4, 0.5, 0.1])
                )
exp = np.arange(2,-17,-1,  dtype=float)
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
    
#%% Plot finite difference test
plt.figure(figsize=(10,4))

# plot CL e CD
plt.subplot(1,2,1)
plt.loglog(h,FDtest_CL, label= r'$C_L$')
plt.loglog(h,FDtest_CD, label = r'$C_D$')

plt.axvline(1e-7,color ='k', linestyle = '--')

plt.xlim(max(h), min(h))
plt.yscale('log')

plt.xlabel('h')
plt.ylabel('FD test')
plt.title('(a)')

plt.grid()
plt.legend()

# Plot residuos
plt.subplot(1,2,2)
plt.loglog(h,FDtest_res)
plt.axvline(1e-7,color ='k', linestyle = '--')

plt.xlim(max(h), min(h))
plt.yscale('log')

plt.xlabel('h')
plt.ylabel('FD test')
plt.title('(b)')


plt.grid()
plt.legend([f'res {i+1}' for i in range(len(res_lltd))], ncols = 2)

plt.tight_layout()

plt.savefig(imagdir / f'FDtest.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
# plt.show()
# %% Dot product test

# Outputs
res_llt_RAD = np.zeros_like(inputs['twist'])
cl_RAD = 0
cd_RAD = 0
# Output seeds
res_lltb = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
clb = 0.7
cdb = 0.8

# Make copies of the output seeds because the RAD
# code will modify them, and this will break the
# dot-product test that we will do later on
res_lltb_copy = res_lltb.copy()
clb_copy = clb
cdb_copy = cdb

# Initialize input seeds for the derivative accumulation
twistb = np.zeros_like(inputs['twist'])
gamab = np.zeros_like(inputs['gama'])

# input seeds
input_seed_RAD = dict(
    res_lltb = res_lltb_copy,
    clb = clb_copy,
    cdb = cdb_copy,
    twistb = twistb,
    gamab = gamab
)

# Call function
llt_b.llt_main_b(**inputs, **input_seed_RAD, res_llt = res_llt_OR, cl = CL_OR,cd = CD_OR)

# Dot product


# %%
