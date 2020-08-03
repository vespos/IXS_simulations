import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
plt.rcParams['image.cmap'] = 'afmhot'
plt.rcParams['figure.figsize'] = (7, 5)

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

import mono_main as mono


bl = mono.build_beamline()
beams = mono.run_process(bl)

# plot zprime histogram
for b in beams:
    beam = beams[b]
    fig, axes = plt.subplots(nrows=1, ncols=2)
    axes[0].set_title(b)
    zprime = raycing.get_zprime(beam)
    axes[0].hist(zprime, bins=100)
    axes[0].set_xlabel('zprime')
    E = beam.E
    axes[1].hist2d(E, zprime, bins=100)
    axes[1].set_xlabel('energy')
    axes[1].set_ylabel('zprime')
    plt.tight_layout()
    saveName = './plots/'+b+'_zprime.png'
    plt.savefig(saveName)
plt.show()