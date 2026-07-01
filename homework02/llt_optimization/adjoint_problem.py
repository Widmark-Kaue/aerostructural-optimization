#%% libs
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import root

from lltlib import LiftingLineOpt

from llt import llt # type: ignore - Código original
from llt_b import llt_b # type: ignore - Código reverso

#%% Inputs
inputs = dict(
            twist = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            span = 8.0,
            chord = 1.0,
            cl0 = 0.0,
            cla = 2*np.pi,
            alpha = 5*np.pi/180,
            vinf = 10,
            rho_air = 1.225,
            )

#%% funcs
def resfun(gama:np.ndarray):
    res,CL,CD = llt.llt_main(gama = gama, **inputs)
    return res

def adjfun(resb:np.ndarray, gama:np.ndarray,func:str):
    # initialize output seeds
    clb = 1
    cdb = 0
    if func.lower() == 'cd':
        clb = 0
        cdb = 1
    
    # Initialize input seeds
    twistb = np.zeros_like(resb)
    gamab = np.zeros_like(resb)

    # Initialize output variables 
    CL = 0
    CD = 0
    res = np.zeros_like(resb)
    
    # call reverse code
    llt_b.llt_main_b(
                    twistb = twistb, 
                    gama = gama, 
                    gamab =gamab,
                    res_llt = res,
                    res_lltb = resb.copy(),
                    cl= CL,
                    clb=clb,
                    cd= CD,
                    cdb = cdb,
                    **inputs)
    
    
    return twistb, gamab

def resAdjfun(resb:np.ndarray, gama:np.ndarray,func:str):
    _,  gamab = adjfun(resb,gama,func)
    adj_res = gamab
    return adj_res
    
#%% Solve physics
gama0  = np.ones_like(inputs['twist'])
sol = root(resfun,gama0)
gama = sol.x
res,CL,CD = llt.llt_main(gama = gama, **inputs)

print('Solve physics:')
print(f'CL = {CL}')
print(f'CD = {CD}')
print(f'gama = {gama}')

# %% Solve adjoint equation for CL
resAdjfun_CL = lambda resb: resAdjfun(resb,gama,'CL')
psi0 = np.ones_like(gama)
sol = root(resAdjfun_CL,psi0)
psiCL = sol.x
print()
print('Adjoint vars CL:')
print(f'psiCL = {psiCL}')

#  Total derivative for CL
dCL_dtwist,_ = adjfun(psiCL,gama,'CL')
print()
print('Total Derivative CL:')
print(f'dCL_dtwist = {dCL_dtwist}')

# %% Solve adjoint equation for CD
resAdjfun_CD = lambda resb: resAdjfun(resb,gama,'CD')
psi0 = np.ones_like(gama)
sol = root(resAdjfun_CD,psi0)
psiCD = sol.x
print()

print('Adjoint vars CD:')
print(f'psiCD = {psiCD}')

#  Total derivative for CD
dCD_dtwist,_ = adjfun(psiCD,gama,'CD')
print()
print('Total Derivative CD:')
print(f'dCD_dtwist = {dCD_dtwist}')

# %% Validate llt class
lltclass = LiftingLineOpt(nvortices=6)
CL,CD, gama = lltclass.solve_llt(inputs['twist']) #type: ignore
dCL_dtwist_class, psiCL_class = lltclass.grad_llt(inputs['twist'],'CL') #type: ignore
dCD_dtwist_class, psiCD_class = lltclass.grad_llt(inputs['twist'],'CD') #type: ignore

print()
print('VALIDATE LLT CLASS:')
print('Solve physics:')
print(f'CL = {CL}')
print(f'CD = {CD}')
print(f'gama = {gama}')


print()
print('Adjoint vars CL:')
print(f'psiCL_class = {psiCL_class}')
print('Total Derivative CL:')
print(f'dCL_dtwist_class = {dCL_dtwist_class}')

print()
print('Adjoint vars CD:')
print(f'psiCD_class = {psiCD_class}')
print('Total Derivative CD:')
print(f'dCD_dtwist_class = {dCD_dtwist_class}')
