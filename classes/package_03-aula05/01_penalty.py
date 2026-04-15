'''
09_penalty.py
uses the quadratic barrier to solve constrained optimization
'''

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# INPUTS

# Define penalty parameter
pen = 0.1

# Define plot bounds
xmin = -1.0
xmax= 1.5
ymin = -0.2
ymax = 1.05

# Define starting point
x0 = -0.5

# Define design functions

def objfun(x):
    f = 1.0 - np.exp(-0.5*(x+0.5)**2)
    return f

def confun(x):
    g1 = -0.3*x

    g2 = 0.3*x-0.3

    return g1,g2

def quad_penalty(x):

    f = objfun(x)
    g1, g2 = confun(x)

    # Compute penalty function
    phi = #???

    return phi

# EXECUTION

results = minimize(quad_penalty, x0=x0, method='Nelder-Mead')
xopt = results.x
print(results)

# PLOTS

fig,ax = plt.subplots(3,1,sharex=True,figsize=(10,6))

# Gather values for the plot
x = np.linspace(xmin, xmax, 1000)
f = []
g1 = []
g2 = []
phi = []

for x_curr in x:
    f_curr = objfun(x_curr)
    g1_curr, g2_curr = confun(x_curr)
    phi_curr = quad_penalty(x_curr)

    f.append(f_curr)
    g1.append(g1_curr)
    g2.append(g2_curr)
    phi.append(phi_curr)

# Plot objective function
ax[0].plot(x,f,'k')

ax[0].plot([0,0],[ymin,ymax],linewidth=0.5,color='gray')
ax[0].plot([1,1],[ymin,ymax],linewidth=0.5,color='gray')

ax[0].set_ylim([ymin, ymax])
ax[0].set_xlim([xmin, xmax])
ax[0].set_xticks([])

ax[0].set_ylabel(r'$f$',fontsize=20)
ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)
ax[0].spines['bottom'].set_visible(False)

# Plot constraints
ax[1].plot(x,g1,'r',x,g2,'r')
ax[1].plot([xmin,xmax],[0,0],linewidth=0.5,color='gray')

ax[1].text(-0.35,0.18,r'$g_1$',fontsize=20,color='r',picker=True)
ax[1].text(1.12,0.15,r'$g_2$',fontsize=20,color='r',picker=True)

ax[1].plot([xmin,xmax],[0,0],linewidth=0.5,color='gray')
ax[1].plot([0,0],[ymin,ymax],linewidth=0.5,color='gray')
ax[1].plot([1,1],[ymin,ymax],linewidth=0.5,color='gray')

ax[1].set_ylim([ymin, ymax*0.5])
ax[1].set_xlim([xmin, xmax])
ax[1].set_xticks([])

ax[1].set_ylabel(r'$g$',fontsize=20)
ax[1].spines['right'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].spines['bottom'].set_visible(False)

# Plot pseudo-objective function
ax[2].plot(x,phi)
ax[2].plot([xopt],[quad_penalty(xopt)],'o')

ax[2].plot([0,0],[ymin,ymax],linewidth=0.5,color='gray')
ax[2].plot([1,1],[ymin,ymax],linewidth=0.5,color='gray')

ax[2].set_ylim([ymin, ymax])
ax[2].set_xlim([xmin, xmax])

ax[2].set_ylabel(r'$\phi$',fontsize=20)
ax[2].set_xlabel(r'$x$',fontsize=20)
ax[2].spines['right'].set_visible(False)
ax[2].spines['top'].set_visible(False)

plt.show()
