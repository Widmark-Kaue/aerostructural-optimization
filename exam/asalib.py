#%%
import numpy as np
from scipy.optimize import root

from copy import deepcopy
from itertools import chain
from typing import Union
from dataclasses import dataclass, field


from asa_module import asa_module as asa # type:ignore
from asa_module_d import asa_module_d as asa_d # type:ignore
from asa_module_b import asa_module_b as asa_b # type:ignore
#%%
@dataclass
class ASAOptimization:
    # Numeric parameters
    npanels: int                                                     # [-]
    pKS: float = 200.0                                               # [-]
    
    # Geometric parameters
    span: float = 8.0                                                # [m]
    chords: Union[np.ndarray, float] = 1.0                           # [m]
    
    # Aerodynamic parameters
    Vinfm: float = 60.0                                              # [m/s]
    rhoAir: float = 1.225                                            # [kg/m^3]
    cd0: float = 0.0270                                              # [-]
    cl0: Union[np.ndarray, float] = 0.0                              # [-]
    cla: Union[np.ndarray, float] = 6.283                            # [1/rad]
    
    # Structural parameters
    E: float = 73.1e9                                                # [Pa]
    rhoMat: float = 2780.0                                           # [kg/m^3]
    sigmaY: float = 324e6                                            # [Pa]
    r: Union[np.ndarray, float] = 0.1                                # [m]
    
    # Operational parameters
    alpha: float = np.deg2rad(5)                                     # [rad]
    fixedMass: float = 700.0                                         # [kg]
    g: float = 9.8                                                   # [m/s^2]
    endurance: float = 4.0 * 60.0 * 60.0                             # [s]
    TSFC: float = 0.5 / 3600.0                                       # [1/s]
    loadFactor: float = 3.0 * 1.5                                    # [-]

    # Utils Dictionary and lists
    outLabel: list = field(init=False, default_factory=list)
    outLabel_b: list = field(init=False, default_factory=list)
    _INPUTS: dict = field(init=False, default_factory=dict)
    _x_hist: list = field(init=False, default_factory=list)
    _f_hist: list = field(init=False, default_factory=list)

    def __post_init__(self):
        self._build_input_dict()
        self._build_outlabel_lists()
        
    #### Properties
    @property
    def x_hist(self):
            return np.array(self._x_hist)
        
    @property
    def f_hist(self):
        return np.array(self._f_hist)

    #### Hidden functions
    
    def _build_input_dict(self):
        # Build vectors
        ones = np.ones(self.npanels, dtype=float)
        chords = self.chords * ones
        cl0 = self.cl0 * ones
        cla = self.cla * ones
        r = self.r * ones
        
        vinf = np.array([self.Vinfm * np.cos(self.alpha), 0, self.Vinfm * np.sin(self.alpha)])
        
        # Geometric dict
        geo = dict(
            span = self.span,
            chords = chords
        )
        
        # Aerodynamic dict
        aero = dict(
            cd0 = self.cd0,                     
            cl0 = cl0, 
            cla = cla,                      
            vinf = vinf,  
            rhoair = self.rhoAir                                                  
        )

        # Structural dict
        struct = dict(
            r = r,                                
            e = self.E,                                                        # [Pa]
            rhomat = self.rhoMat,                                              # [kg/m^3]
            sigmay = self.sigmaY,                                              # [Pa]
        )
        
        # Operational parameters
        oper = dict(
            fixedmass = self.fixedMass,                                        # [kg]
            g = self.g,                                                        # [m/s^2]
            endurance = self.endurance,                                        # [s]
            tsfc = self.TSFC,                                                  # [1/s]
            loadfactor = self.loadFactor
        )
        # Numeric dict
        num = dict(pks = self.pKS)
        
        self._INPUTS = geo | aero | struct | oper | num
        
    def _build_outlabel_lists(self):
        self.outLabel = ['resllt', 'resfem', 'liftexcess', 'margins', 'ksmargin', 'fb', 'weight', 'sref', 'cl']
        # names for output seeds Reverse AD
        self.outLabel_b = [f'{self.outLabel[i]}b' for i in range(7)]

    def _run_asa(self, twist: np.ndarray, gama: np.ndarray, t: np.ndarray, d: np.ndarray) -> dict:
        """
        Runs the original code (asa_main) and returns the output variables mapped to a dictionary.
        
        Parameters:
        -----------
        gama : np.ndarray
            Circulation distribution.
        twist : np.ndarray
            Twist distribution.
        t : np.ndarray
            Thickness distribution.
        d : np.ndarray
            Displacements vector.
            
        Returns:
        --------
        dict
            Dictionary with the outputs mapped to outLabel keys:
            ['resllt', 'resfem', 'liftexcess', 'margins', 'ksmargin', 'fb', 'weight', 'sref', 'cl']
        """
        # Create dictionary of the design and state variables
        dv = {
            'gama': gama,
            'twist': twist,
            't': t,
            'd': d
        }
        # Merge design variables with the pre-built parameter dictionary
        inputs = self._INPUTS | dv
        
        # Call the original Fortran solver
        outputs = asa.asa_main(**inputs)
        
        # Map the outputs to their corresponding labels
        return dict(zip(self.outLabel, outputs))

    def _run_asa_b(self, twist: np.ndarray, gama: np.ndarray, t: np.ndarray, d: np.ndarray, 
                  output_seeds: dict) -> dict:
        """
        Runs the reverse AD code (asa_main_b) to compute gradients.
        
        Parameters:
        -----------
        gama : np.ndarray
            Circulation distribution.
        twist : np.ndarray
            Twist distribution.
        t : np.ndarray
            Thickness distribution.
        d : np.ndarray
            Displacements vector.
        output_seeds : dict
            Dictionary of the seeds for the outputs. Keys should be in outLabel_b,
            e.g., {'reslltb': ..., 'weightb': ..., etc.}.
            Any missing seeds will be automatically initialized to 0.
            
        Returns:
        --------
        dict
            Dictionary containing the computed input gradients:
            {'gamab': np.ndarray, 'twistb': np.ndarray, 'tb': np.ndarray, 'db': np.ndarray}
        """
        # 1. Initialize original outputs as zero
        outputs_reverse = dict(
            resllt = np.zeros(self.npanels, dtype=float),
            resfem = np.zeros(2*(self.npanels+1), dtype = float),
            liftexcess = 0.0,
            margins = np.zeros(2*self.npanels, dtype=float),
            ksmargin = 0.0,
            fb = 0.0,
            weight = 0.0,
            sref = 0.0,
            cl = 0.0
        )
        
        
        # 2. Initialize output seeds
        out_b = deepcopy(output_seeds)
                
        # 3. Initialize input seeds (accumulators for reverse AD)
        inp_b = dict(
            gamab= np.zeros_like(gama, dtype=float),
            twistb= np.zeros_like(twist, dtype=float),
            tb= np.zeros_like(t, dtype=float),
            db= np.zeros_like(d, dtype=float)
        )
        
        # 4. Prepare design and state variables
        dv = dict(
            gama= gama,
            twist= twist,
            t= t,
            d= d
        )
        
        # 5. Build inputs for reverse AD call
        inputs_b = deepcopy(self._INPUTS) | dv | inp_b | outputs_reverse | out_b
        
        # Call the Fortran reverse sweep solver
        asa_b.asa_main_b(**inputs_b)
        
        # Return the computed input gradients
        return inp_b
        
    def _resfunc(self, desVars:np.ndarray,stateVars:np.ndarray):
        #Get design vars
        twist = desVars[:self.npanels]
        t = desVars[self.npanels:]
        
        #Get state vars
        gama = stateVars[:self.npanels]/0.1
        d = stateVars[self.npanels:]/100.0
        
        out = self._run_asa(twist=twist,gama=gama,t =t, d =d)
        res = np.hstack([out['resllt'], out['resfem']])
        return res
    
    def _adjfunc(self, desVars:np.ndarray, stateVars:np.ndarray,resb:np.ndarray, func:str):
        #Get design vars
        twist = desVars[:self.npanels]
        t = desVars[self.npanels:]
        
        #Get state vars
        gama = stateVars[:self.npanels]
        d = stateVars[self.npanels:]
        
        #Initialize output seeds
        reslltb = resb[:self.npanels]
        resfemb = resb[self.npanels:]
        out_b = [reslltb,resfemb] + [0]*5
        output_seeds = dict(zip(self.outLabel_b, out_b))
        output_seeds['marginsb'] = np.zeros(2*self.npanels,dtype=float)
        output_seeds[f'{func.lower()}b'] = 1 # activate interest function
        
        # Call reverse code (Compute input seeds)
        input_seeds = self._run_asa_b(twist=twist,gama=gama,t=t,d=d, output_seeds=output_seeds)
        
        stateVarsb = np.hstack([input_seeds['gamab'], input_seeds['db']])
        desVarsb = np.hstack([input_seeds['twistb'], input_seeds['tb']])
        return stateVarsb, desVarsb
            
    def _resAdjfunc(self,desVars:np.ndarray, stateVars:np.ndarray, resb:np.ndarray, func:str):
        stateVarsb,_ = self._adjfunc(desVars,stateVars, resb,func)
        adj_res = stateVarsb
        return adj_res
    
    #### Public functions
    
    def clean_history(self):
        self._x_hist = []
        self._f_hist = []
    
    def solve_asa(self,  desVars:np.ndarray):
        # Design variables
        desVars[self.npanels:] =desVars[self.npanels:]/100
        twist = desVars[:self.npanels]
        t = desVars[self.npanels:]
        
        # Solve physics
        resfunc = lambda stateVars: self._resfunc(desVars, stateVars)
        # Use professor's default initial guesses (gama = 80.0, d = 0.001)
        gama0 = np.array([80.0, 80.0, 80.0, 80.0])
        d0 =  np.array([0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001])
        # Scale the initial guess for the solver's normalized variables
        stateVars0 = np.hstack([gama0 * 0.1, d0 * 100.0])
        sol = root(resfunc, stateVars0, options={'xtol': 1e-6})
        
        # Split state vars (applying same scaling: gama = stateVars/0.1, d = stateVars/100)
        gama = sol.x[:self.npanels]/0.1
        d = sol.x[self.npanels:]/100.0
        
        # Call original code
        out = self._run_asa(twist=twist,gama= gama, t = t, d =d)
        
        return out
    
       
        
        
        