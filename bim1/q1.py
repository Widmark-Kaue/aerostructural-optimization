#%% Imports
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from utils import set_aiaa_style

#%% plot settings

set_aiaa_style(16)

#%% Image settings
imagdir = Path('.')
if not imagdir.name == 'bim1':
    imagdir = Path('.','bim1', 'images')

imagdir.mkdir(exist_ok=True)

dpi = 600
format = 'pdf'
saveflag = True
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
