import numpy as np
from scipy.optimize import root
from dataclasses import dataclass, field

from llt import llt # type: ignore - Código original
from llt_d import llt_d # type: ignore - Código direto
from llt_b import llt_b # type: ignore - Código reverso

@dataclass
class LiftingLineOpt:
    CL_target:float = 0.5         # [-]
    span:float = 8.0              # [m]
    chord:float = 1.0             # [m]             
    cl0: float = 0.0              # 
    cla: float = 2*np.pi          # [1/rad]
    alpha: float = 5*np.pi/180    # [rad]
    Vinf: float  =  10            # [m/s]
    rho_air: float = 1.225        # [kg/m^3]             
    save_history:bool = False
    _x_hist:list = field(init=False, default_factory=list)
    _f_hist:list = field(init=False, default_factory=list)
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
        
    @property
    def x_hist(self):
        return np.array(self._x_hist)
    
    @property
    def f_hist(self):
        return np.array(self._f_hist)

    def _clean_history(self):
        self._x_hist = []
        self._f_hist = []
    
    def _resfun(self, twist:np.ndarray,gama:np.ndarray):
        res, CL, CD = llt.llt_main(twist=twist,gama =gama, **self.inputs)
        return res
    
    def _adjfun(self, twist:np.ndarray,gama:np.ndarray,resb:np.ndarray,func:str):
        # Initialize output seeds
        clb = 1
        cdb = 0
        if func.lower() == 'cd':
            clb = 0
            cdb = 1
        
        # Initialize input seeds
        twistb = np.zeros_like(twist)
        gamab = np.zeros_like(gama)

        # Initialize output variables 
        CL = 0
        CD = 0
        res = np.zeros_like(resb)
        
        # call reverse code
        llt_b.llt_main_b(
                        twist = twist, 
                        twistb = twistb, 
                        gama = gama,
                        gamab = gamab,
                        cl = CL, 
                        clb =clb, 
                        cd = CD, 
                        cdb = cdb,
                        res_llt = res, 
                        res_lltb = resb.copy(),
                         **self.inputs)
        
        
        return twistb, gamab
    
    def _resAdjfun(self, twist:np.ndarray,gama:np.ndarray,resb:np.ndarray,func:str):
        _,  gamab = self._adjfun(twist,gama,resb,func)
        adj_res = gamab
        return adj_res

    
    def solve_llt(self,twist:np.ndarray):
        resfun = lambda gama: self._resfun(twist,gama)
        gama0  = np.zeros_like(twist)
        
        # Solve llt
        sol = root(resfun,gama0)  
        gama = sol.x
        
        # Compute coefficients
        _,CL,CD = llt.llt_main(twist=twist,gama=gama,**self.inputs)
        
        return CL, CD, gama
          
    def grad_llt(self, twist:np.ndarray,func:str):
        # Solve physics
        _,_,gama = self.solve_llt(twist)
        
        # Solve adjoint equation
        resAdjfun = lambda resb: self._resAdjfun(twist,gama,resb,func)
        psi0 = np.ones_like(twist)
        sol = root(resAdjfun,psi0)
        psi = sol.x
        
        # Compute total derivative
        dfunc_dtwist,_ = self._adjfun(twist,gama,psi,func)
        
        return dfunc_dtwist, psi
        
    
    def objfun(self,twist:np.ndarray):
        # Solve physics
        _,CD,_ = self.solve_llt(twist)
        
        if self.save_history:
            self._x_hist.append(twist)
            self._f_hist.append(CD)
        return CD

    def objfungrad(self,twist:np.ndarray):
        dCD_dtwist,_ = self.grad_llt(twist,'CD')
        return dCD_dtwist
    
    
    def confun(self, twist:np.ndarray):
        # Solve physics
        CL,_,_ = self.solve_llt(twist)
        h = CL - self.CL_target
        return h
    
    def confungrad(self, twist:np.ndarray):
        dCL_dtwist,_ = self.grad_llt(twist,'CL')
        return dCL_dtwist