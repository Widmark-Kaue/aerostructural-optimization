#%% libs
import numpy as np

from scipy.optimize import root
from asalib import ASAOptimization
#%% Create asa opt object
asaObj = ASAOptimization(4)

#%% Values to validate
Gama_val = np.array([12.99803608,14.2879062,14.2879062,12.99803608])
d_val = np.array([ 2.69493469e-02,-8.96491634e-03,9.57299675e-03,-7.85795127e-03,
            -9.88815861e-45,-4.23707119e-67,9.57299675e-03,7.85795127e-03,
            2.69493469e-02,8.96491634e-03])

psi_A_val = np.array([-92.62945248, -16.9979671, -16.8732846, -115.59593901])
psi_S_val = np.array([1.76989719e-04, -4.41232719e-05, 8.87470713e-05, -4.41174285e-05,
                      -1.53738245e+02, -8.13436714e-08, 8.87470714e-05, 4.41174285e-05,
                      1.76989720e-04, 4.41232719e-05])

dKSmargin_dtwist_val = np.array([-2.63697857, -1.16166158, -1.16166158, -2.63697857])
dKSmargin_dt_val = np.array([1.11410770e-03, 6.61735632e+01, 6.61735633e+01, 1.11410770e-03])

#%% 5.4 Solve physics
twist = np.array([0.0,0.0,0.0,0.0])
t = np.array([0.005,0.005,0.005,0.005])
desVars = np.hstack([twist, t])

resfunc = lambda stateVars: asaObj._resfunc(desVars=desVars,stateVars = stateVars)

gama0 = np.array([80.0,80.0,80.0,80.0])
d0 =  np.array([0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001])

stateVars0 = np.hstack([gama0, d0])
print(stateVars0)
sol = root(resfunc,  stateVars0, options={'xtol': 1e-6})
# sol = root(resfunc,  stateVars0)
print('SOLVE  PHYSICS:')
print(sol)

gama = sol.x[:asaObj.npanels]/0.1
d = sol.x[asaObj.npanels:]/100.0
stateVars = np.hstack([gama,  d])
# gama
match = np.allclose(gama, Gama_val)
print(f'Gama: Code =  {gama}')
print(f'Gama: 5.4 Values = {Gama_val}')
print(f'Match Gama: {match}\n')
print()
match = np.allclose(d, d_val)
print(f'd: Code =  {d}')
print(f'd: 5.4 Values = {d_val}')
print(f'Match d: {match}\n')

#%% 5.5 Solve Adjoint problem
reslltb = np.ones_like(twist)*1e2
resfemb = np.ones_like(d)
resb0 = np.hstack([reslltb, resfemb])

resAdjfunc = lambda resb: asaObj._resAdjfunc(desVars, stateVars,resb,'ksmargin')
# resAdjfunc(resb0)
print('SOLVE ADJOINT PROBLEM')
sol = root(resAdjfunc, resb0)
print(sol)

psi_ksmargin = sol.x
psi_A = psi_ksmargin[:asaObj.npanels]
psi_S = psi_ksmargin[asaObj.npanels:]
print('Adjoint vars:')
match_psi_A = np.allclose(psi_A, psi_A_val)
print(f'psi_A: Code =  {psi_A}')
print(f'psi_A: 5.5 Values = {psi_A_val}')
print(f'Match psi_A: {match_psi_A}\n')

match_psi_S = np.allclose(psi_S, psi_S_val)
print(f'psi_S: Code =  {psi_S}')
print(f'psi_S: 5.5 Values = {psi_S_val}')
print(f'Match psi_S: {match_psi_S}\n')

_,dKSmargin_dx = asaObj._adjfunc(desVars,stateVars,psi_ksmargin, 'ksmargin')
dKSmargin_dtwist = dKSmargin_dx[:asaObj.npanels]
dKSmargin_dt = dKSmargin_dx[asaObj.npanels:]

print('Total derivatives:')
match_dtwist = np.allclose(dKSmargin_dtwist, dKSmargin_dtwist_val)
print(f'dKSmargin_dtwist: Code =  {dKSmargin_dtwist}')
print(f'dKSmargin_dtwist: 5.5 Values = {dKSmargin_dtwist_val}')
print(f'Match dKSmargin_dtwist: {match_dtwist}\n')

match_dt = np.allclose(dKSmargin_dt, dKSmargin_dt_val)
print(f'dKSmargin_dt: Code =  {dKSmargin_dt}')
print(f'dKSmargin_dt: 5.5 Values = {dKSmargin_dt_val}')
print(f'Match dKSmargin_dt: {match_dt}\n')

# %%
