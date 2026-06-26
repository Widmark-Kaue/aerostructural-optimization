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
    rhoKS: float = 10   
    F: np.ndarray = field(init=False, default_factory=lambda: np.array([]))
    a: np.ndarray = field(init=False,default_factory=lambda: np.array([]))
    
    def __post_init__(self):
        
        # build force vectors
        self._build_force_vector() 

        # build weight matrix
        self._build_weight_matrix()

    @property
    def alpha(self):
        return self.R*self.E/self.sigma_y
    
    def _build_force_vector(self):
        self.F = np.array([self.F1, 0, self.F2, 0]).reshape(-1,1)
    
    def _build_weight_matrix(self):
        self.a = np.array([ [-6,  2,  0,  0],
                            [-6,  4,  0,  0],
                            [6,   4,  -6, 2],
                            [6,   2,  -6, 4]])
        
    
    def _build_stiffness_matrix(self,ta,tb):
        # Construindo Matriz de rigidez
        cte = np.pi*self.R**3*self.E/1e3
        mat = np.diag([12*ta+12*tb, 4*ta+4*tb, 12*tb, 4*tb])

        mat[0,1:] = np.array([-6*ta+6*tb,  -12*tb, 6*tb])
        mat[1,2:] = np.array([-6*tb, 2*tb])
        mat[2,3] = -6*tb

        K = cte*(mat + mat.T - np.eye(mat.shape[0])*mat)
        return K

    def _resfun(self, A,b,c):
        return A@b - c 

    def _find_state_vars(self,x):
        K = self._build_stiffness_matrix(x[0],x[1])
        resfun = lambda p: self._resfun(K,p,self.F)
        
        p0 = np.ones_like(self.F)
        sol = root(resfun,p0)
        p = sol.x
        return p
    
    def _find_adjoint_vars(self,delrdelp,delfdelp):
        resfun = lambda psi: self._resfun(delrdelp.T,psi,-delfdelp.T)

        psi0= np.ones_like(delfdelp)
        sol = root(resfun,psi0)
        psi = sol.x
        return psi
        
    def objfun(self, x):
        ta,tb = x
        m = 2*np.pi*self.R*self.rho/1e3 *(ta + tb)
        return m
    
    def objfungrad(self, x):
        cte =  2*np.pi*self.R*self.rho/1e3
        dmdx = cte*np.ones(2)
        return dmdx

    def confun(self, x):
        p = self._find_state_vars(x)
        gis = 1 - self.alpha**2 * self.a@p
        return gis
    
    def confunKS(self,x):
        gis = self.confun(x)
        Gks = -1/self.rhoKS * np.log(np.sum(np.exp(-self.rhoKS*gis)))
        return Gks 
    
    def confunKSgrad(self, x):
        p = self._find_state_vars(x)
        si = self.a@p
        delsidelp = self.a
        delGksdelp = -2*self.alpha**2 * (si @ delsidelp)
        return delGksdelp 
        
        
