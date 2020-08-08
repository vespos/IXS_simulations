import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')

from importlib import reload
import bl_test as test
reload(test)

bl_fun = test.build_beamline
process_fun = test.run_process

class align_crystal():
    def __init__(self, bl_fun, process_fun):
        self.bl = bl_fun()
        self.beams = process_fun(self.bl)
        self.n_crystal = np.sum([1 for oe in self.bl.oes if 'crystal' in oe.name])


ca = align_crystal(bl_fun, process_fun)
bl = ca.bl
beam_s = ca.beams[0]

# fig, ax = plt.subplots()
# for ii, beam in enumerate(ca.beams):
#     if beam==0:
#         continue
#     ax.plot(beam.Jsp, '.', label=ii, markersize=1)
# plt.legend()
# plt.show()