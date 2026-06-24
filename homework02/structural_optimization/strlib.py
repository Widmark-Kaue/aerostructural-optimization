import numpy as np
from scipy.optimize import root


def resfun(A,p,F):
    return A@p - F 

def objfun(ta,tb):
    pass
