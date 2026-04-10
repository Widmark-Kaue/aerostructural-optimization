""""
    Question 1 from Homerwork 1: 
        - Plotting the z-component of velocity induced by a horseshoe vortex
"""
#%% Imports
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from utils import set_aiaa_style

#%% plot settings

set_aiaa_style(16)

dpi = 600
format = 'pdf'
saveflag = True
#%% Image Path settings
imagdir = Path('.')
if not imagdir.absolute().name == 'bim1':
    imagdir = Path('.','bim1')

imagdir = imagdir.joinpath('images')
imagdir.mkdir(exist_ok=True, parents=True)
print("Images will be saved in:", imagdir)

#%% Plotting z-component of velocity induced by a horseshoe vortex
L = 1.0  # Length of the domain
Gamma = 1.0  # Circulation

y = np.linspace(-2, 2)
w = Gamma*L/ (4*np.pi) /(y**2-L**2/4)

plt.figure(figsize=(8, 6))
plt.plot(y, w, 'darkblue', linewidth=2)
plt.xlabel('y [m]')
plt.ylabel('w [m/s]')

plt.xlim(-2, 2)
plt.ylim(np.min(w), np.max(w))

plt.grid()
plt.savefig(imagdir / f'q1.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
plt.show() 
