import numpy as np
from scipy.optimize import minimize 
import matplotlib.pyplot as plt

x_hist = []
f_hist = []

def objfun(X):

    x1 = X[0]
    x2 = X[1]
    
    ff = 100*(x2-x1**2)**2 + (np.sqrt(2)-x1)**2
    
    x_hist.append(X)
    f_hist.append(ff)
    
    return ff
    
def objfun_grad(X):

   x1 = X[0]
   x2 = X[1]
   
   dfdx1 = -400*(x2-x1**2)*x1 - 2*(np.sqrt(2)-x1)
   dfdx2 = 200*(x2-x1**2)
   
   grad = np.array([dfdx1, dfdx2])

   return grad

#Chute Inicial
X0 = [0.0, 0.0]

#Otimizador
result = minimize(objfun, X0, method='bfgs', tol=1e-9, jac=objfun_grad)

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
