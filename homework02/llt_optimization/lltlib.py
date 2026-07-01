import numpy as np
from scipy.optimize import root
from dataclasses import dataclass, field

from llt import llt # type: ignore - Código original
from llt_d import llt_d # type: ignore - Código direto
from llt_b import llt_b # type: ignore - Código reverso

@dataclass
class LiftingLineOpt:
    span:float = 8.0              # [m]
    chord:float = 1.0             # [m]             
    cl0: float = 0.0              # 
    cla: float = 2*np.pi          # [1/rad]
    alpha: float = 5*np.pi/180    # [rad]
    Vinf: float  =  10            # [m/s]
    rho_air: float = 1.225        # [kg/m^3]             
    save_history:bool = False
    x_hist:list = field(init=False, default_factory=list)
    f_hist:list = field(init=False, default_factory=list)
    inputs:dict = field(init=False,default_factory=dict)
    
    
    def __post_init__(self):
        
        # build inputs dict
        self._build_inputs_dict() 

    
    def _build_inputs_dict(self):
        self.inputs = dict(
            span = self.span,
            chord = self.chord,
            cl0 = self.cl0,
            cla = self.cla,
            alpha = self.alpha,
            vinf = self.Vinf,
            rho_air = self.rho_air,
            )
        

    def _resfun(self, twist:np.ndarray,gama:np.ndarray):
        res, CL, CD = llt.llt_main(twist=twist,gama =gama, **self.inputs)
        return res

    
    def solve_llt(self,twist:np.ndarray):
        resfun = lambda gama: self._resfun(twist,gama)
        gama0  = np.zeros_like(twist)
        
        # Solve llt
        sol = root(resfun,gama0)  
        
        gama =sol.x
        _,CL,CD = llt.llt_main(twist=twist,gama=gama,**self.inputs)
        
        return CL, CD, gama
          
    def grad_llt(self, twist:np.ndarray,func:str):
        
        CL,CD,gama = self.solve_llt(twist)
        
        res = np.zeros_like(gama)
        
        # Output seeds
        resb = np.zeros_like(gama)
        clb = 1
        cdb = 0
        if func.lower == 'cd':
            clb = 0
            cdb = 1
        
        # Initialize input seeds for derivative accumulation
        twistb = np.zeros_like(twist)
        gamab = np.zeros_like(gama)
        
        # Call function
        llt_b.llt_main_b(
                        twist = twist, twistb = twistb, gama = gama,gamab = gamab,
                         **self.inputs,
                        cl = CL, clb =clb, cd = CD, cdb = cdb,res_llt = res, res_lltb = resb)
        
    # def objfun(self, x):
    #     ta,tb = x
    #     m = 2*np.pi*self.R*self.rho/1e3 *(ta + tb)
    #     if self.save_history:
    #         self.x_hist.append(x)
    #         self.f_hist.append(m)
    #     return m
    
    # def objfungrad(self, x):
    #     cte =  2*np.pi*self.R*self.rho/1e3
    #     dmdx = cte*np.ones(2)
    #     return dmdx

    # def confun(self, x):
    #     # Solve physics
    #     K = self._build_stiffness_matrix(x[0],x[1])
    #     p = self._find_state_vars(K)
        
    #     # Compute isolated constraints
    #     si = self.A@p
    #     gis = 1 - self.alpha**2 * si**2
    #     return gis
    
    # def confungrad(self, x):
    #     # Solve physics
    #     K = self._build_stiffness_matrix(x[0],x[1])
    #     p = self._find_state_vars(K)
        
    #     # Compute si vector
    #     si = self.A@p
        
    #     # Compute partial derivatives
    #     delsidelp = self.A
    #     delrdelp = K
    #     delrdelx = self._build_delrdelx(p)
        
    #     dgidx = []
    #     for i in range(len(si)):
    #         # Compute partial derivate of g related to state vars (for each constraint)
    #         delgidelp = -2*self.alpha**2*si[i]*delsidelp[i]
            
    #         # Solve adjoint equation for each constraint
    #         psigis = self._find_adjoint_vars(delrdelp,delgidelp)
           
    #         # Total derivative for each constraint
    #         dgidx.append(psigis.T @ delrdelx)
        
    #     return dgidx