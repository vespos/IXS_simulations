import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')

import spectro_9keV_v1 as spec9
import spectro_45keV as spec4

bl9 = spec9.build_beamline()
bl4 = spec4.build_beamline()

print('9kev')
spec9.beamline_stats(bl9)
print('\n')
print('4.5kev')
spec9.beamline_stats(bl4)

""" propagation direction """
angles = [0, spec9.mirr1_out, spec9.xtal1_out, spec9.xtal2_out, spec9.mirr2_out]
n = []
for angle in angles:
    n.append(np.array([0,
        np.cos(angle),
        np.sin(angle)])) # unitary vector along the beam direction (global coord)

bls = [bl4, bl9]
coord = []
for bl in bls:
    coord = np.asarray([
        bl.source.center[1:],
        bl.mirr1.center[1:],
        bl.xtal1.center[1:],
        bl.xtal2.center[1:],
        bl.mirr2.center[1:],
        bl.s4.center[1:]
    ])

    plt.plot(coord[:,0], coord[:,1], '-o', markersize=10)
    labels = ['S','M1','C1','C2','M2','det']
    for lab, c in zip(labels,coord):
        plt.annotate(lab, (c[0]+10, c[1]+10))
plt.xlim(-50,1100)
plt.ylim(-200,80)
plt.savefig('./layout_2.png')
plt.show()