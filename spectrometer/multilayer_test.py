import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sin, cos, tan, arctan, degrees, radians

sys.path.append(r'C:\Users\espov\Documents\Python\xrt')
from xrt.backends.raycing.materials import Material
from xrt.backends.raycing.materials import Multilayer
from elements import LMultilayer
# ---------------------------------------------------------------------------

angle_i = np.arange(0.1,3,0.001) # grazing angle of incidence
beamInDotNormals = sin(radians(angle_i))

# MLs and amplitudes
E0 = 8798
mSi = Material('Si', rho=2.33)
mRu = Material('Ru', rho=11.0)
mC = Material('C', rho=2.26)
# p = 0.35
# ml = LMultilayer(tLayer=mC, tThickness=30, bLayer=mRu, bThickness=30,
#     nPairs=150, substrate=mSi, p=p, E0=E0)
ml1 = Multilayer(tLayer=mC, tThickness=13.15, bLayer=mRu, bThickness=13.15,
    nPairs=150, substrate=mSi)
A1 = np.asarray([ml1.get_amplitude(E0, ang)[0] for ang in beamInDotNormals])

E0 = 4399
mSi = Material('Si', rho=2.33)
mNi = Material('Ni', rho=8.9)
mC = Material('C', rho=2.26)
ml2 = Multilayer(tLayer=mC, tThickness=30, bLayer=mNi, bThickness=30,
    nPairs=100, substrate=mSi)
A2 = np.asarray([ml2.get_amplitude(E0, ang)[0] for ang in beamInDotNormals])

# plots
plt.figure()
plt.title('ML reflectivity')
plt.plot(angle_i, np.abs(A1)**2, label='9keV')
plt.plot(angle_i, np.abs(A2)**2, label='4.5keV')
plt.xlabel('Incident angle (deg)')
plt.ylabel('Reflectivity')
plt.yscale('log')
plt.legend()
plt.show()