import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')

import spectro_9keV as spec

bl = spec.build_beamline()
beams = spec.run_process(bl)

plt.hist(bl.source.make_energy(), bins=100)
plt.show()