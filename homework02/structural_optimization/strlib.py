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
    save_history:bool = False
    F: np.ndarray = field(init=False, default_factory=lambda: np.array([]))
    a: np.ndarray = field(init=False,default_factory=lambda: np.array([]))
    x_hist:list = field(init=False, default_factory=list)
    f_hist:list = field(init=False, default_factory=list)
    
    def __post_init__(self):
        
        # build force vectors
        self._build_force_vector() 

        # build weight matrix
        self._build_weight_matrix()

    @property
    def alpha(self):
        return self.R*self.E/self.sigma_y
    
    def _build_force_vector(self):
        self.F = np.array([self.F1, 0, self.F2, 0])
    
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
    
    def _build_delrdelx(self,p):
        cte = np.pi*self.R**3*self.E/1000
        v1,phi1,v2,phi2 = p
        mat = np.array([[12*v1-6*phi1, 12*v1+6*phi1-12*v2+6*phi2    ],
                        [-6*v1+4*phi1, 6*v1+4*phi1-6*v2+2*phi2      ],
                        [0,            -12*v1-6*phi1+12*v2-6*phi2   ],
                        [0,             6*v1+2*phi1-6*v2+4*phi2     ]
                        ])
        delrdelx = cte*mat
        return delrdelx

    def _resfun(self, A,b,c):
        return A@b - c 

    def _find_state_vars(self,K):
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
        if self.save_history:
            self.x_hist.append(x)
            self.f_hist.append(m)
        return m
    
    def objfungrad(self, x):
        cte =  2*np.pi*self.R*self.rho/1e3
        dmdx = cte*np.ones(2)
        return dmdx

    def confun(self, x):
        K = self._build_stiffness_matrix(x[0],x[1])
        p = self._find_state_vars(K)
        si = self.a@p
        gis = 1 - self.alpha**2 * si**2
        return gis
    
    def confungrad(self, x):
        K = self._build_stiffness_matrix(x[0],x[1])
        p = self._find_state_vars(K)
        si = self.a@p
        delsidelp = self.a
        delrdelp = K
        delrdelx = self._build_delrdelx(p)
        
        dgidx = []
        for i in range(len(si)):
            delgidelp = -2*self.alpha**2*si[i]*delsidelp[i]
            psigis = self._find_adjoint_vars(delrdelp,delgidelp)
            dgidx.append(psigis.T @ delrdelx)
        
        return dgidx
            
    def confunKS(self,x):
        gis = self.confun(x)
        Gks = -1/self.rhoKS * np.log(np.sum(np.exp(-self.rhoKS*gis)))
        return Gks 
    
    def confunKSgrad(self, x):
        K = self._build_stiffness_matrix(x[0],x[1])
        p = self._find_state_vars(K)
        si = self.a@p
        delsidelp = self.a
        delGksdelp = -2*self.alpha**2 * (si @ delsidelp)
        delrdelp = K
        psiG = self._find_adjoint_vars(delrdelp,delGksdelp)
        
        delrdelx = self._build_delrdelx(p)
        dGksdx = psiG.T @ delrdelx
        return dGksdx
        
        
