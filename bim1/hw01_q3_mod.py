# IMPORTS
import numpy as np

#==========================

def rosenbrock(x):
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
    
    return

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

    return

#==========================

def nm_opt(n):
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

    return

#==========================

def de_opt(n):
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

    return

#==========================

def cg_opt(n):
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

    return

#==========================

def bfgs_opt(n):
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

    return

#==========================
