import numpy as np
from scipy.optimize import minimize, differential_evolution
import matplotlib.pyplot as plt

x_hist = []
f_hist = []

def objfun(X):

    x1 = X[0]
    x2 = X[1]
    
    ff = 36 + x1**2 + x2**2 - 8*x1 -10*np.cos(2*np.pi*x1) - 10*np.cos(2*np.pi*x2)
    
    x_hist.append(X)
    f_hist.append(ff)
    
    return ff
    
bounds = [[-5,5], [-5,5]]

result = differential_evolution(objfun, bounds, seed=1, maxiter=5000, polish=False)

print(result)

xopt= result.x

print('xopt:', xopt)

fig = plt.figure()
plt.plot(x_hist)
plt.xlabel('nfev')
plt.ylabel('X')

fig = plt.figure()
plt.semilogy(f_hist)
plt.xlabel('nfev')
plt.ylabel('f')

plt.show()
