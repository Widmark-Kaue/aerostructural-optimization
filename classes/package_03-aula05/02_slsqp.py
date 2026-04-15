'''
02_slsqp.py
Constrained Optimization
Cap. Ney Sêcco 2020
'''

# IMPORTS
import numpy as np
from scipy.optimize import minimize
import auxmod as am
import matplotlib.pyplot as plt

# EXECUTION
xlist = []
flist = []
g1list = []
g2list = []

# Define objective function

def objfun(x):
    f = 1-np.exp(-0.1*(x[0]**2 + x[1]**2))
    flist.append(f)
    xlist.append(x.copy())
    return f

def objfungrad(x):
    gradf = np.array([np.exp(-0.1*(x[0]**2 + x[1]**2))*0.2*x[0],
                      np.exp(-0.1*(x[0]**2 + x[1]**2))*0.2*x[1]])
    return gradf

# Define constraints

def confun1(x):
    g = x[0]+x[1]-2
    g1list.append(g)

    return g

def confun1grad(x):
    gradg = np.array([1,
                      1])
    return gradg

def confun2(x):
    g = x[1]-1
    g2list.append(g)
    
    return g

def confun2grad(x):
    gradg = np.array([0,
                      1])
    return gradg

# Create list of constraints 
con1 = { 'type': 'ineq', 
        'fun': confun1,
        'jac': confun1grad}

con2 = { 'type': 'ineq', 
        'fun': confun2,
        'jac': confun2grad}

cons = [con1, con2]

# Define starting point
x0 = np.array([3.0, 2.0])

# Define optimization bounds
bounds = [(-4,4),(-4,4)]

# Run optimizer
result = minimize(objfun, x0, jac = objfungrad, bounds=bounds, constraints=cons, method='SLSQP')

# Print results

print(result)                       
xopt = result.x

"""
    Nota-se que a segunda restrição é inativa visto na opção multipliers o segundo parâmetro é 0.
    Se retirarmos a segunda restrição da lista, o mínimo continua igual. 
"""


# Plot optimization history
fig = plt.figure()
plt.subplot(311)
plt.plot(xlist,'o-')
plt.ylabel('x',fontsize=20)
plt.subplot(312)
plt.plot(flist,'o-')
plt.ylabel('f',fontsize=20)
plt.subplot(313)
plt.plot(g1list,'o-')
plt.plot(g2list,'o-')
plt.plot([0,len(g1list)-1],[0,0],'gray',linewidth=0.5)
plt.ylabel('g',fontsize=20)
plt.xlabel('evaluations',fontsize=20)
plt.tight_layout()

# Save path before it is modified
xk = xlist.copy()

# Plot results
fig = plt.figure()
ax = plt.gca()
am.plot_contour(objfun, ax,
                xmin=-3*1.3, xmax=3*1.3, ymin=-3, ymax=3, zmin=10**-3, zmax=10**(-0.1))
am.plot_contour(confun1, ax,
                xmin=-3*1.3, xmax=3*1.3, ymin=-3, ymax=3, zmin=0, zmax=0, nlevels=None) # type: ignore
am.plot_contour(confun2, ax,
                xmin=-3*1.3, xmax=3*1.3, ymin=-3, ymax=3, zmin=0, zmax=0, nlevels=None) # type: ignore
am.plot_path(ax, xk, xopt=xopt)

plt.show()

