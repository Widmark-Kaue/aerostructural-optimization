'''
ks_test.py
This function demonstrates how KS function works
'''

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt

# EXECUTION

# Define constraint functions
sigma = 1.0
g1 = lambda x: -x+2
g2 = lambda x: (x - 2)*(x - 8)
g3 = lambda x: 5*2.0/(np.sqrt(3*sigma)*np.pi**0.33)*(1.0 - ((x-5)/sigma)**2)*np.exp(-((x-5)/sigma)**2)-3.0

# Define sampling points
n = 100
xmin = 0.0
xmax = 10.0
x = np.linspace(xmin, xmax, n)

# Compute function at sampling points
y1 = g1(x)
y2 = g2(x)
y3 = g3(x)

# Define KS function parameter
rhoKS = 10.0

# Compute KS function
gks = 1/rhoKS * np.log(np.exp(rhoKS*y1) + np.exp(rhoKS*y2) + np.exp(rhoKS*y3)) # CHANGE THIS LINE

# Plot results
fig = plt.figure()
plt.plot(x, y1, 'k')
plt.plot(x, y2, 'k')
plt.plot(x, y3, 'k')
plt.plot(x, np.zeros(n),'gray',linewidth=0.5)
plt.plot(x, gks, 'r')
plt.xlabel('x',fontsize=15)
plt.ylabel('g',fontsize=15)
plt.xlim([xmin, xmax])
plt.show()
