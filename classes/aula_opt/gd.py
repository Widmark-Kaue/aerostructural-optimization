import numpy as np
from scipy.optimize import minimize, differential_evolution
import matplotlib.pyplot as plt

x_hist = []
f_hist = []

def objfun(X):

    x1 = X[0]
    x2 = X[1]
    
    ff = x1**2 - 2*x1*x2 + np.sqrt(2)*x2**2
    
    x_hist.append(X)
    f_hist.append(ff)
    
    return ff
    
def objfun_grad(X):
 
   x1 = X[0]
   x2 = X[1]
   
   dfdx1 = 2*x1 - 2*x2
   dfdx2 = -2*x1 + 2*np.sqrt(2)*x2
   
   grad = np.array([dfdx1, dfdx2])
   
   return grad

X0 = np.array([2,0])

maxiter = 300

X = X0.copy()

for ii in range(maxiter):
     
     ff = objfun(X)
     
     grad = objfun_grad(X)
     
     if np.linalg.norm(grad) < 1e-6:
     	break
     alpha = 0.02
     
     pk = -grad
     
     X =  X +alpha*pk



xopt= X

print('xopt:', xopt)

fig = plt.figure()
plt.plot(x_hist)
plt.xlabel('nfev')
plt.ylabel('X')

fig = plt.figure()
plt.plot(f_hist)
plt.xlabel('nfev')
plt.ylabel('f')

plt.show()
