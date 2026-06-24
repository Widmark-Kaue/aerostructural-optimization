import numpy as np
from scipy.optimize import root
from dataclasses import dataclass, field

@dataclass
class StructuralOpt:
    R:float             # [m]
    E:float             # [Pa]             
    sigma_y: float      # [Pa]
    rho: float          # [kg/m^3]
    F1: float           # [N]
    F2: float           # [N]
    rhoKs: float = 10   
    F: np.ndarray = field(init=False, default_factory=lambda: np.array([]))
    A: np.ndarray = field(init=False, default_factory=lambda: np.array([]))
    
    def __post_init__(self):
        # build force vectors
        self.F = np.array([self.F1, 0, self.F2, 0])
         

    def _build_stiffness_matrix(self):
        pass

    def _resfun(self, p):
        return self.A@p - self.F 

    def objfun(self, ta,tb):
        pass
