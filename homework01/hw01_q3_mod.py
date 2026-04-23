# IMPORTS
import numpy as np
import os
from scipy.optimize import minimize, differential_evolution
#==========================

def rosenbrock(x:np.ndarray) -> float:
    '''
    Computes the multidimension Rosenbrock function
    at point x. The number of dimensions is equal to
    the length of x.
    
    INPUTS
    x : array(n) -> Point where the function will be evaluated.
    must be a numpy array.
    
    OUTPUTS
    1) float : Rosenbrock function value
    
    EXAMPLE
    rosenbrock(np.array([1.0, 1.0])) -> 0.0
    rosenbrock(np.array([-1.0, 0.0, 0.3])) -> 114.0
    '''
    
    # ADD CODE HERE
    n =  len(x)
    f = np.zeros(n-1)
    
    for i in range(n-1):
        f[i] = 100*(x[i+1] - x[i]**2)**2 + (1-x[i])**2
    
    return np.sum(f)

#==========================

def rosenbrock_grad(x):
    '''
    Computes the gradient of multidimension Rosenbrock function
    at point x. The number of dimensions (n) is equal to
    the length of x.
    
    INPUTS
    x : array(n) -> Point where the function gradient will be evaluated.
    must be a numpy array.
    
    OUTPUTS
    1) array(n) : Rosenbrock function gradient
    
    EXAMPLE
    rosenbrock_grad(np.array([1.0, 1.0])) -> np.array([0.0, 0.0])
    rosenbrock_grad(np.array([-1.0, 0.0, 0.3])) -> np.array([-404., -202.,   60.])
    '''
    
    # ADD CODE HERE
    n = len(x)
    grad = np.zeros(n)
    
    grad[0] = -400*x[0]*(x[1] - x[0]**2) - 2*(1-x[0])
    grad[n-1] = 200*(x[n-1] - x[n-2]**2)
    
    for j in range(1, n-1):
        grad[j] = 200*(x[j] - x[j-1]**2) - 400*x[j]*(x[j+1] - x[j]**2) - 2*(1-x[j]) 
    
    return grad

#==========================

def nm_opt(n, tol = 1e-12):
    '''
    Optimizes the function using the Nelder-Mead Simplex Algorithm,
    starting at the point x = [0, 0, ..., 0].

    INPUTS
    n : integer -> Number of dimensions of the Rosenbrock function.

    OUTPUTS
    1) array : Optimum values of design variables (xopt).
    2) float : Optimum value of the objective function (fopt)
    3) integer : Number of function evaluations (nfev)

    EXAMPLES:
    nm_opt(2) -> np.array([1.00000028, 1.00000051]), 3.5147904406946316e-13, 169
    nm_opt(3) -> np.array([0.99999983, 0.99999961, 0.99999923]), 4.3261171377916544e-13, 338
    '''

    # ADD CODE HERE
    x = np.zeros(n)
    options = {'maxiter': 200000, 'fatol': tol, 'adaptive': True}
    results = minimize(rosenbrock, x, method='Nelder-Mead', options=options)

    x_opt = results.x
    f_opt = results.fun
    nfev = results.nfev
    
    print(f'nm_opt: {results.message}')
    
    return x_opt, f_opt, nfev

#==========================

def de_opt(n, tol = 1e-12):
    '''
    Optimizes the function using the Differential Evolution Algorithm,
    starting at the point x = [0, 0, ..., 0].

    INPUTS
    n : integer -> Number of dimensions of the Rosenbrock function.

    OUTPUTS
    1) array : Optimum values of design variables (xopt).
    2) float : Optimum value of the objective function (fopt)
    3) integer : Number of function evaluations (nfev)

    EXAMPLES:
    de_opt(2) -> np.array([0.99999997, 0.99999995]), 6.686425291204817e-15, 1710
    de_opt(3) -> np.array([0.99999989, 0.9999998 , 0.99999959]), 9.727695757591168e-14, 4680
    '''

    # ADD CODE HERE
    bounds = [(-4, 4)] * n
    results = differential_evolution(
        rosenbrock, 
        bounds, 
        seed = 1,  #type: ignore
        maxiter=5000, 
        polish = False, 
        atol=tol) #type: ignore
    x_opt = results.x
    f_opt = results.fun
    nfev = results.nfev
    
    print(f'de_opt: {results.message}')
    return x_opt, f_opt, nfev

#==========================

def cg_opt(n, tol = 1e-6):
    '''
    Optimizes the function using the Conjugate Gradient Algorithm,
    starting at the point x = [0, 0, ..., 0].

    INPUTS
    n : integer -> Number of dimensions of the Rosenbrock function.

    OUTPUTS
    1) array : Optimum values of design variables (xopt).
    2) float : Optimum value of the objective function (fopt)
    3) integer : Number of function and gradient evaluations (nfev+njev)

    EXAMPLES:
    cg_opt(2) -> np.array([0.99999949, 0.99999897]), 2.639508773268233e-13, 108
    cg_opt(3) -> np.array([1., 1., 1.]), 8.796853489285454e-21, 134
    '''

    # ADD CODE HERE
    x = np.zeros(n)
    options = {'maxiter': 200000}
    results = minimize(rosenbrock, x, method='CG', jac=rosenbrock_grad, tol = tol,  options=options)
    x_opt = results.x
    f_opt = results.fun
    nfev = results.nfev
    njac = results.njev
    
    print(f'cg_opt: {results.message}')
    return x_opt, f_opt, nfev + njac

#==========================

def bfgs_opt(n, tol = 1e-6):
    '''
    Optimizes the function using the BFGS Algorithm,
    starting at the point x = [0, 0, ..., 0].

    INPUTS
    n : integer -> Number of dimensions of the Rosenbrock function.

    OUTPUTS
    1) array : Optimum values of design variables (xopt).
    2) float : Optimum value of the objective function (fopt)
    3) integer : Number of function and gradient evaluations (nfev+njev)

    EXAMPLES:
    bfgs_opt(2) -> np.array([1., 1.]), 2.3804051301997937e-18, 50
    bfgs_opt(3) -> np.array([1., 1., 1.00000001]), 3.8670228678944245e-17, 76
    '''

    # ADD CODE HERE
    x = np.zeros(n)
    options = {'maxiter': 200000}
    results = minimize(rosenbrock, x, method='BFGS', jac=rosenbrock_grad, tol = tol, options=options)
    x_opt = results.x
    f_opt = results.fun
    nfev = results.nfev
    njac = results.njev
    
    print(f'bfgs_opt: {results.message}')
    return x_opt, f_opt, nfev + njac

#==========================

def test_func(func, inputs, expected_outputs):

    # Run the desired function
    current_outputs = func(*inputs)

    # Compare output
    # test_result = (current_outputs== expected_outputs)

    # Print log
    print('')
    print('======================')
    print('Testing function ' + func.__name__)
    print('  Inputs:')
    print(inputs)
    print('  Current outputs:')
    print(current_outputs)
    print('  Expected outputs:')
    print(expected_outputs)
    # print(' Test successful?')
    # print(test_result)
    # print('======================')
    print('')

    # Return success flag
    # return test_result


if __name__ == "__main__":
    # Test the functions
    os.system('clear') # for linux
    print(50 * '-')
    print('Testing rosenbrock functions:')
    test_func(rosenbrock, [np.array([1.0, 1.0])], 0.0)
    test_func(rosenbrock, [np.array([-1.0, 0.0, 0.3])], 114.0)
    test_func(rosenbrock_grad, [np.array([1.0, 1.0])], np.array([0.0, 0.0]))
    test_func(rosenbrock_grad, [np.array([-1.0, 0.0, 0.3])], np.array([-404., -202.,   60.]))
    

    print('\n\n'+50 * '-')
    print('Testing optimization functions:')
    
    test_func(nm_opt, [2], (np.array([1.00000028, 1.00000051]), 3.5147904406946316e-13, 169))
    test_func(nm_opt, [3], (np.array([0.99999983, 0.99999961, 0.99999923]), 4.3261171377916544e-13, 338))

    test_func(de_opt, [2], (np.array([0.99999997, 0.99999995]), 6.686425291204817e-15, 1710))
    test_func(de_opt, [3], (np.array([0.99999989, 0.9999998 , 0.99999959]), 9.727695757591168e-14, 4680))
        
    test_func(cg_opt, [2], (np.array([0.99999949, 0.99999897]), 2.639508773268233e-13, 108))
    test_func(cg_opt, [3], (np.array([1., 1., 1.]), 8.796853489285454e-21, 134))
    
    test_func(bfgs_opt, [2], (np.array([1., 1.]), 2.3804051301997937e-18, 50))
    test_func(bfgs_opt, [3], (np.array([1., 1., 1.00000001]), 3.8670228678944245e-17, 76))
    
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    
    plt.figure()
    x1 = np.linspace(-2, 2, 100)
    x2 = np.linspace(-1, 3, 100)
    X1, X2 = np.meshgrid(x1, x2)
    comb = np.array([X1.flatten(), X2.flatten()]).T
    Z = np.zeros_like(X1).flatten()
    for i in range(len(X1.flatten())):
        Z[i] = rosenbrock(comb[i])
    Z = Z.reshape(X1.shape)
    plt.contourf(X1, X2, Z, levels = 150, cmap = 'viridis', norm = LogNorm())
    
    plt.xlabel('$X_1$')
    plt.ylabel('$X_2$')
    
    plt.colorbar(label = 'Rosenbrock function value')
    plt.show() 