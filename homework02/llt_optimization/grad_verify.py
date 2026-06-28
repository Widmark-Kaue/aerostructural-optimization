#%% libs
import numpy as np
import matplotlib.pyplot as plt

import copy

from llt import llt # type: ignore - Código original
from llt_d import llt_d # type: ignore - Código direto
from llt_b import llt_b # type: ignore - Código reverso
#%% Inputs
inputs = dict(
            twist = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            gama = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
            span = 8.0,
            chord = 1.0,
            cl0 = 0.0,
            cla = 2*np.pi,
            alpha = 5*np.pi/180,
            vinf = 10,
            rho_air = 1.225,
            )


inputs_seed = dict(
                twistd = np.array([0.3, 0.1, 0.2, 0.4, 0.6, 0.5]),
                gamad = np.array([0.2, 0.3, 0.6, 0.4, 0.5, 0.1])
                )

#%% Finite Diference test
h = 10**(np.arange(-1,-16,-1))
inputs1 = copy.copy(inputs)

inputs1['twist'] = inputs['twist']*(1+h)
inputs1['gama'] = inputs['gama']*(1+h)

res_llt,CL,CD= llt.llt_main(**inputs)
res_llt1,CL1,CD1= llt.llt_main(**inputs1)

res_lltd_FD = (res_llt1-res_llt)/h
print(res_llt)
print(res_llt1)
print(res_lltd_FD)

# %%


# %%
