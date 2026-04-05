import numpy as np
import matplotlib.pyplot as plt

L = 1.0  # Length of the domain
Gamma = 1.0  # Circulation

y = np.linspace(-2, 2)
w = Gamma*L/ (4*np.pi) /(y**2-L**2/4)

plt.figure()
plt.plot(y, w)
plt.xlabel('y')
plt.ylabel('w')

plt.show() 
