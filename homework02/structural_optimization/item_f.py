import numpy as np
from scipy.optimize import root
from strlib import StructuralOpt

# Constantes
E = 73.1e9          # [Pa]
sigma_y = 324e6     # [Pa]
rho = 2780          # [kg/m^3]
R = 0.05            # [m]
F1 = F2 = 1000      # [N]
rhoKs = 10          # 
ta = tb = 5         # [mm]

# Construindo Matriz de rigidez
cte = np.pi*R**3*E/1e3
mat = np.diag([12*ta+12*tb, 4*ta+4*tb, 12*tb, 4*tb])

mat[0,1:] = np.array([-6*ta+6*tb,  -12*tb, 6*tb])
mat[1,2:] = np.array([-6*tb, 2*tb])
mat[2,3] = -6*tb

K = cte*(mat + mat.T - np.eye(mat.shape[0])*mat)

print(f'K = {K}')

# Vetor de forças
F = np.array([F1, 0, F2, 0])

# Resolver residuo
r = lambda d: K @ d - F
d0  = np.ones_like(F)

res = root(r,d0)
p = res.x
# print(f'res = {res}')
print(f'p = {p}')

# Derivadas parciais
v1, phi1, v2, phi2 = p

drdx = cte*np.array([[12*v1 - 6*phi1, 12*v1 + 6*phi1 - 12*v2 + 6*phi2   ],
                 [-6*v1 + 4*phi1, 6*v1 + 4*phi1 - 6*v2 + 2*phi2     ],
                 [0             , -12*v1 - 6*phi1 + 12*v2 - 6*phi2  ],
                 [0             , 6*v1 + 2*phi1 - 6*v2 + 4*phi2     ]
                 ])
print(f'drdx = {drdx}')

# Calculate s_i and g_i
si = np.array([
    -6*v1 + 2*phi1,
    -6*v1 + 4*phi1,
    6*v1 + 4*phi1 - 6*v2 + 2*phi2,
    6*v1 + 2*phi1 - 6*v2 + 4*phi2
])

alpha = R * E / sigma_y
gis = 1 - alpha**2 * si**2

# Calculate weights
w = np.exp(-rhoKs * gis)
wi = w / np.sum(w)

# ds_i/dp
delsidelp = np.array([
    [-6,  2,  0,  0],
    [-6,  4,  0,  0],
    [ 6,  4, -6,  2],
    [ 6,  2, -6,  4]
])

# dGksdp
dGksdp = -2 * alpha**2 * (wi * si) @ delsidelp
print(f'dGksdp = {dGksdp}')

# Derivada total m
dmdx = 2*np.pi*R*rho/1000 * np.ones(2)
print(f'dmdx {dmdx}')

# Derivada total Gks
rg = lambda psi: K.T @ psi + dGksdp.T # residuo da eq adjunta para Gks
psi0 = np.ones((4,1))
res = root(rg, psi0)

psiG = res.x
print(f'psiG = {psiG}')
dGksdx = psiG.T @ drdx
print(f'dGksdx = {dGksdx}')


#%% Validação da classe
prob = StructuralOpt(R, E,sigma_y,rho,F1,F2)
x=[ta,tb]

print('---------------------------')
print('Validação da classe:')
print('dmdx= ', prob.objfungrad(x))
print('dGksdx = ',  prob.confunKSgrad(x))
