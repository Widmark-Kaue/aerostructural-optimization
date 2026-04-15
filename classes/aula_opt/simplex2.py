import numpy as np
from scipy.optimize import minimize 
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
    
#Chute Inicial
X0 = [3.9, 0.0]

#Otimizador
result = minimize(objfun, X0, method='Nelder-Mead', tol=1e-14)

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
