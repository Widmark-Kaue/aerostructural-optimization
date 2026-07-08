#%%
import numpy as np
from typing import Optional,Union

from dataclasses import dataclass, field

from asa_module import asa_module as asa # type:ignore
from asa_module_d import asa_module_d as asa_d # type:ignore
from asa_module_b import asa_module_b as asa_b # type:ignore
#%%
@dataclass
class ASAOptimization:
    # Numeric
    npanels:int
    pKS:float = 200
    
    # Geometric parameters
    span: float = 8 # [m]
    chords: Union[np.ndarray, float] = 1.0 # [m]
    
    # Aerodynamic parameters
    Vinfm: float = 60 # [m/s]
    rhoAir: float  = 1.225 # [kg/m^3]
    cd0: float = 0.0270
    cl0: Union[np.ndarray, float] = 0
    cla: Union[np.ndarray, float] = 6.283
    
    # Structural parameters
    E: float = 73.1e9 # [Pa]
    rhoMat:float = 2780.0 # [kg/m^3]
    sigmaY:float = 324e6 # [Pa]
    
    r:Union[np.ndarray, float] = 0.1 # [m]
    
    # Operational parameters
    alpha:float = np.rad2deg(5) # [rad]
    fixedMass: float  = 700.0 # [kg]
    g: float  = 9.8 # [m/s^2]
    endurance = 4.0*60.0*60.0,                                      # [s]
    TSFC = 0.5/3600.0,                                              # [1/s]
    loadFactor = 3.0*1.5

    # Utils Dictionary and lists
    _INPUTS:dict = field(init=False, default_factory=dict)
    _outlabel: list = field(init=False,default_factory=list)
    _outlabel_b: list = field(init=False, default_factory=list)
    _x_hist:list = field(init=False, default_factory=list)
    _f_hist:list = field(init=False, default_factory=list)

    def __post_init__(self):
        self._build_input_dict()
        self._build_outlabel_lists()
    
    def _build_input_dict(self):
        # Build vectors
        ones = np.ones(self.npanels, dtype=float)
        chords = self.chords*ones
        cl0 = self.cl0*ones
        cla = self.cla*ones
        r = self.r*ones
        
        vinf = np.array([self.Vinfm*np.cos(self.alpha), 0, self.Vinfm*np.sin(self.alpha)])
        
        # Geometric dict
        geo = dict(
               span = self.span,
               chords = self.chords
            )
        
        # Aerodynamic dict
        aero =  dict(
            cd0 = self.cd0,                     
            cl0 = cl0, 
            cla = cla,                      
            vinf= vinf,  
            rhoair = self.rhoAir                                                  
        )

        # Structural dict
        struct = dict(
            r = r,                                
            e = self.E,                                                       #[Pa]
            rhomat = self.rhoMat,                                                #[kg/m^3]
            sigmay = self.sigmaY,                                                 #[Pa]
            )
        
        # Operational parameters
        oper = dict(
            fixedmass = self.fixedMass,                                              # [kg]
            g = self.g,                                                        # [m/s2]
            endurance = self.endurance,                                      # [s]
            tsfc = self.TSFC,                                              # [1/s]
            loadfactor = self.loadFactor
        )
        # Numeric dict
        num = dict(pks = self.pKS)
        
        self._INPUTS = geo | aero | struct | oper
        
    def _build_outlabel_lists(self):
        self.outLabel = ['resllt','resfem','liftexcess','margins','ksmargin','fb','weight','sref','cl']
        # names for output seeds Reverse AD
        self.outLabel_b = [f'{self.outLabel[i]}b' for i in range(7)]