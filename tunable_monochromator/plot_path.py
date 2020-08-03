import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
plt.rcParams['image.cmap'] = 'afmhot'
plt.rcParams['figure.figsize'] = (7, 5)

from importlib import reload

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')

import xrt.backends.raycing as raycing
import xrt.plotter as xrtplot
import xrt.runner as xrtrun

from xrt.plotter import XYCAxis, XYCPlot
from xrt.backends.raycing.sources import GeometricSource
from xrt.backends.raycing.screens import Screen
from xrt.backends.raycing.oes import OE
from xrt.backends.raycing.materials import Crystal, CrystalSi

import mono_9keV_Si111 as mono
reload(mono)

bl = mono.build_beamline()
beams = mono.run_process(bl)
det_beam = beams['detector']

filtx = np.logical_and(det_beam.x>-0.02, det_beam.x<0.02)
filtz = np.logical_and(det_beam.z>-0.07, det_beam.z<-0.03)
filt = np.logical_and(filtx,filtz)
if filt.sum()==0:
    filt = np.ones(det_beam.E.shape, dtype=bool)
    print("\nNo filter\n")

en = det_beam.E
opt_path = det_beam.path
fit = np.polyfit(en, opt_path, 1)
fit = np.poly1d(fit)
yfit = fit(en)

plt.figure()
plt.title('Optical path full detector')
plt.hist2d(en, opt_path-yfit, bins=150)
# plt.hist2d(en, opt_path, bins=200)
plt.plot(en, np.zeros(en.shape), linewidth=1.5, color='black')
plt.xlabel('Energy (eV)')
plt.ylabel('Optical path (mm)')
plt.tight_layout()

en_filt = det_beam.E[filt]
opt_path_filt = det_beam.path[filt]
fit = np.polyfit(en_filt, opt_path_filt, 1)
fit = np.poly1d(fit)
yfit = fit(en_filt)

plt.figure()
plt.title('Optical path ROI')
plt.hist2d(en_filt, opt_path_filt-yfit, bins=150)
# plt.hist2d(en, opt_path, bins=200)
plt.plot(en, np.zeros(en.shape), linewidth=1.5, color='black')
plt.xlabel('Energy (eV)')
plt.ylabel('Optical path (mm)')
plt.tight_layout()

plt.figure()
plt.title('Energy dependence x')
plt.hist2d(en, det_beam.x, bins=200)
plt.xlabel('Energy (eV)')
plt.ylabel('x (mm)')
plt.tight_layout()

plt.figure()
plt.title('Energy dependence z')
plt.hist2d(en, det_beam.z, bins=200)
plt.xlabel('Energy (eV)')
plt.ylabel('z (mm)')
plt.tight_layout()

plt.show()