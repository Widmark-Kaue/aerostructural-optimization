import numpy as np
from scipy.optimize import root

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

cte_G = 2*(R*E/sigma_y)**2
dGksdp = cte_G*np.array([
                    4*36*v1 - 2*36*v2 + 36*phi2,
                    40*phi1 - 36*v2 + 16*phi2,
                    -72*v1 - 36*phi1 + 72*v2 - 36*phi2,
                    36*v1 + 16*phi1 - 36*v2 + 20*phi2
                    ])
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

