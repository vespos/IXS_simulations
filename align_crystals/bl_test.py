import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
# plt.switch_backend('agg')

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')
import os

import xrt.backends.raycing as raycing
import xrt.plotter as xrtplot
import xrt.runner as xrtrun

from xrt.plotter import XYCAxis, XYCPlot
from xrt.backends.raycing.sources import GeometricSource
from xrt.backends.raycing.screens import Screen
from xrt.backends.raycing.oes import OE
import xrt.backends.raycing.materials as rm

E0 = 14000
dE = 0.1
d_source = 200
FWHM_x = 1000
FWHM_z = 300
source_kwargs = {
    'name': 'source',
    'center': [0,-d_source,0],
    'nrays': 1e5,
    'dx': FWHM_x/1000/2.634, # mm
    'dz': FWHM_z/1000/2.634, # mm
    'dxprime': 0.03/1e6, # rad
    'dzprime': 0.03 /1e6, # rad
    'distE': 'flat', # lines, normal, flat
    'energies': [E0-dE, E0+dE],
    'polarization': 'h'
    # 'polarization': 'v'
    }
angle = 0

Si_kwargs = {
    'hkl': [3,3,3],
    # 'tK': 300 # [K]
    'tK': 120 # [K]
    }
Si_crystal = rm.CrystalSi(**Si_kwargs)
Si_theta = Si_crystal.get_Bragg_angle(E0)

n = np.array([0, np.cos(angle), np.sin(angle)])
center = [0,0,0]
pitch = Si_theta
crystal0_kwargs = {
    'name': 'crystal0',
    'center': center,
    'pitch': 'auto',
    'positionRoll': 0,
    'material': Si_crystal,
    'alpha': 0
    }
angle = angle + (2*Si_theta - 2*Si_crystal.get_dtheta(E0, 0))

center = crystal0_kwargs['center'] + np.array([0, np.cos(angle), np.sin(angle)])*200
crystal1_kwargs = {
    'name': 'crystal1',
    'center': center,
    'pitch': 'auto',
    'positionRoll': np.pi,
    'material': Si_crystal,
    'alpha': 0,
    }


def build_beamline():
    bl = raycing.BeamLine()

    # Source and optical elements
    bl.source = GeometricSource(bl=bl, **source_kwargs)
    oes = []
    oes.append(OE(bl=bl, **crystal0_kwargs))
    oes.append(OE(bl=bl, **crystal1_kwargs))
    # bl.oes = oes

    # Screens
    # print('\nScreens position:')
    # for kwargs in screens_kwargs:
    #     Screen(bl=bl, **kwargs)
    #     print('\t' + kwargs['name'] + ': ' + str(bl.screens[-1].center))
    # print('\n')
    return bl

def run_process(bl):
    beam = [bl.source.shine()]
    outDict = {}
    for ii, oe in enumerate(bl.oes):
        # beam_s = screen.expose(beam=beam[ii])
        if 'slits' in oe.name:
            beam.append(oe.propagate(beam=beam[ii], needNewGlobal=True)[0])
            # beam.append(oe.propagate(beam=beam[ii][0], needNewGlobal=True))
        elif 'crl' in oe.name:
            beam.append(oe.multiple_refract(beam=beam[ii])[0])
            # beam.append(oe.multiple_refract(beam=beam[ii][0]))
        else:
            beam.append(oe.reflect(beam=beam[ii])[0])
            # beam.append(oe.reflect(beam=beam[ii][0]))
        # outDict[screens.name] = beam_s
            # note: two beam are output by the functions 'reflect'/'propagate'/...: 
            # [0] is the beam in global coordinates and [1] in local coordinates.
            # The global beam should be passed to subsequent element.
            # The local output can be suppressed by kwarg 'needLocal=False'
    # screen = bl.screens[-1]
    # beam_s = screen.expose(beam=beam[-1])
    # outDict[screen.name] = beam_s
    return beam
