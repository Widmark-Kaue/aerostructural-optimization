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

#%% 5.4 Solve physics
twist = np.array([0.0,0.0,0.0,0.0])
t = np.array([0.005,0.005,0.005,0.005])
designVars = np.hstack([twist, t])

resfunc = lambda stateVars: asaObj._resfunc(designVars=designVars,stateVars = stateVars)

gama0 = np.array([80.0,80.0,80.0,80.0])
d0 =  np.array([0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001])

stateVars0 = np.hstack([gama0,d0])
print(stateVars0)
sol = root(resfunc,  stateVars0)
print(sol)

gama = sol.x[:asaObj.npanels]/100
d = sol.x[asaObj.npanels:]/0.1
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
resfemb = np.ones_like(t)
resb0 = np.hstack([reslltb, resfemb])

resAdjfunc = lambda resb: asaObj._resAdjfunc(designVars, stateVars,resb,'ksmargin')


