'''
vtol.py
Test case for hand-implementation of adjoint method for a VTOL aircraft
'''

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root

#=======================================
# USER-INPUTS

# Acceleration of gravity
g = 9.81

# Air density
rho = 1.1

# Desired speed
V = 12

# Aircraft weight
W = 4.84*g

# Wing area
S = 0.7

# Drag polar (CD = CD0 + K*CL^2)
CD0 = 0.0251
K = 0.1324

#=======================================
# DEFINITIONS

# Dynamic pressure
q = 0.5*rho*V**2

def main(p,x):

    # Break design and state variables
    sigma = x[0]
    T = p[0]
    CL = p[1]

    # Compute residuals
    r1 = T*np.sin(sigma) + q*S*CL - W
    r2 = T*np.cos(sigma) - q*S*(CD0 + K*CL**2)
    res = [r1, r2]

    # Compute design functions
    f = T

    return res, f

# Define function to return residuals only
def residuals(p,x):

    res, f = main(p,x)

    return res

#=======================================
# Adjoint solution for gradient-based optimization

def objfun(x):

    # Solve the direct problem to get state variables
    p0 = [W, 0.5] # Guess for state variables
    result = root(residuals, p0, args=(x))
    p = result.x

    # Get corresponding objective function
    res, f = main(p,x)

    return f

def objfun_grad(x):

    # Solve the direct problem to get state variables
    p0 = [W, 0.5] # Guess for state variables
    result = root(residuals, p0, args=(x))
    p = result.x

    # Break input variables
    sigma = x[0]

    # Break state variables
    T = p[0]
    CL = p[1]

    # Compute the terms of the adjoint equation
    dfdx = np.array([0])

    dfdp = np.array([1, 0])

    drdp = np.array([[np.sin(sigma), q*S],
                     [np.cos(sigma), -2*q*S*K*CL]])

    drdx = np.array([T*np.cos(sigma), -T*np.sin(sigma)])

    # Solve the adjoint equation
    psi = np.linalg.solve(drdp.T, - dfdp)

    # Compute gradient
    DfDx = dfdx + psi.T @ drdx

    return DfDx

# Run optimization
x0 = [0.0]
result = minimize(objfun, x0, jac=objfun_grad, method='bfgs')
sigmaotm = result.x[0]

# Guess for state variables
p0 = [W, 0.5]

# Solve for state variables
p0 = [W, 0.5] # Guess for state variables
result = root(residuals, p0, args=(result.x))
p = result.x
Totm = p[0]
CLotm = p[1]

print('Adjoint optimization')
print('CLotm:',CLotm)
print('sigmaotm:',sigmaotm*180/np.pi)
print('Totm:',Totm)
print('')

#=======================================
# PLOTS

# Plot bounds
sigmamin = 0.0
sigmamax = 45.0
Tmin = 6.0
Tmax = 8.0
zmin = Tmin
zmax = Tmax
nlevels = 10

# GENERATE OBJECTIVE FUNCTION CONTOUR
delta = 100
x = np.linspace(sigmamin, sigmamax, delta)
y = np.linspace(Tmin, Tmax, delta)
X, Y = np.meshgrid(x, y)
nx,ny = X.shape
Z = Y # Since objfun = T
fig = plt.figure(figsize=(8,6))
ax = fig.gca()
plt.contour(X, Y, Z, np.linspace(zmin, zmax, nlevels))

# Decoration
plt.xlabel(r'$\sigma$ [deg]',fontsize=25)
plt.ylabel(r'$T$ [N]',fontsize=25)
ax = plt.gca()
ax.tick_params(axis='both', labelsize=15)
plt.tight_layout()

# GENERATE LINE OF PHYSICAL VALUES
delta = 100
sigma_list = np.linspace(sigmamin, sigmamax, delta)
T_list = []
for sigma in sigma_list:

    # Guess for state variables
    T0 = 0.3*W
    CL0 = 0.5
    p0 = [T0, CL0]

    # Solve for state variables
    result = root(residuals, p0, args=([sigma*np.pi/180]))
    p = result.x
    T = p[0]
    CL = p[1]

    T_list.append(T)

# Plot physical solutions
plt.plot(sigma_list, T_list, 'r', linewidth=2)
plt.text(sigmamin+0.75*(sigmamax-sigmamin),6.95,r'$\vec{r}=0$',color='r', fontsize=25)

# Plot optimum point
plt.plot(sigmaotm*180/np.pi, Totm, 'ok', markersize=10)
plt.text(sigmaotm*180/np.pi,Totm-0.1,'optimum',color='k', fontsize=25, ha='center', va='center')
plt.savefig('vtol.png')
