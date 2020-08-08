import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
# plt.switch_backend('agg')

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')
import os

import re
from scipy.spatial.transform import Rotation

import xrt.backends.raycing as raycing
import xrt.plotter as xrtplot
import xrt.runner as xrtrun

from xrt.plotter import XYCAxis, XYCPlot
from xrt.backends.raycing.sources import GeometricSource
from xrt.backends.raycing.screens import Screen
from xrt.backends.raycing.oes import OE, ParabolicalMirrorParam, DoubleParaboloidLens
from xrt.backends.raycing.materials import Material, Crystal, CrystalSi
from xrt.backends.raycing.apertures import RectangularAperture, DoubleSlit

import kwargs_mono as kwa

repeat = 10
E0 = 11215
dE = 0.05
# hkl = [6,6,0]
# hkl = [4,4,4]
# hkl = [4,4,0]
hkl = [3,3,3]
withLens = True

# Save folders
subfolder = '{:d}keV/Si{}{}{}/'.format(int(np.round(E0/1000)), hkl[0], hkl[1], hkl[2])
# save_prefix = 'slits_open_'
save_prefix = 'slits_7um_'
save_suffix = ''
do_plot = 'all'


if E0==17795 and hkl==[6,6,0]:
    miscut = np.radians(27.98)
    manual_pitch_corr = [0, 1.65e-5, 0, 0] # 18 keV, 660
elif E0==16000 and hkl==[6,6,0]:
    miscut = np.radians(27.98)
    manual_pitch_corr = [0, 1.06e-5, 0, 0] # 16 keV, 660
elif E0==14000 and hkl==[6,6,0]:
    miscut = np.radians(27.98)
    manual_pitch_corr = [0, 7.9e-6, 0, 0] # 14 keV, 660
elif E0==12000 and hkl==[6,6,0]:
    miscut = np.radians(27.98)
    manual_pitch_corr = [0, 6.4e-6, 0, 0] # 12 keV, 660

elif E0==17795 and hkl==[4,4,4]:
    miscut = np.radians(9.3)
    manual_pitch_corr = [0, 2.8e-6, 0, 0] # 18 keV, 444
elif E0==16000 and hkl==[4,4,4]:
    miscut = np.radians(9.3)
    manual_pitch_corr = [0, 2.7e-6, 0, 0] # 16 keV, 444
elif E0==14000 and hkl==[4,4,4]:
    miscut = np.radians(9.3)
    manual_pitch_corr = [0, 2.6e-6, 0, 0] # 14 keV, 444
elif E0==12000 and hkl==[4,4,4]:
    miscut = np.radians(9.3)
    manual_pitch_corr = [0, 2.6e-6, 0, 0] # 12 keV, 444

elif E0==11215 and hkl==[4,4,0]:
    miscut = np.radians(27.98)
    manual_pitch_corr = [0, 2.875e-5, 0, 0] # 11.215 keV, 440
elif E0==11215 and hkl==[3,3,3]:
    miscut = np.radians(9.3)
    manual_pitch_corr = [0, 4.9e-6, 0, 0] # 11.215 keV, 440
elif E0==11215 and hkl==[4,4,4]:
    miscut = np.radians(9.3)
    manual_pitch_corr = [0, 1.15e-5, 0, 0] # 11.215 keV, 440

else:
    miscut = 0
    manual_pitch_corr = [0, 0, 0, 0]

oes_kwargs, screens_kwargs = kwa.mono_kwargs(E0, dE, hkl, miscut, \
    crystal_pitch_corr=manual_pitch_corr, withLens=withLens)




def build_beamline():
    bl = raycing.BeamLine()

    # Source and optical elements
    bl.source = GeometricSource(bl=bl, **oes_kwargs[0])
    oes = []
    print('OEs position:')
    for kwargs in oes_kwargs[1:]:
        if 'crystal' in kwargs['name']:
            oes.append(OE(bl=bl, **kwargs))
        elif 'parabolic_mirror' in kwargs['name']:
            oes.append(ParabolicalMirrorParam(bl=bl, **kwargs))
        elif 'flat_mirror' in kwargs['name']:
            oes.append(OE(bl=bl, **kwargs))
        elif 'crl' in kwargs['name']:
            oes.append(DoubleParaboloidLens(bl=bl, **kwargs))
        elif 'slits' in kwargs['name']:
            oes.append(RectangularAperture(bl=bl, **kwargs))
        else:
            print('OE not implemented (build beamline)')
        print('\t' + oes[-1].name + ': ' + str(oes[-1].center))
    bl.oes = oes

    # Screens
    print('\nScreens position:')
    for kwargs in screens_kwargs:
        Screen(bl=bl, **kwargs)
        print('\t' + kwargs['name'] + ': ' + str(bl.screens[-1].center))
    print('\n')
    return bl


def run_process(bl):
    beam = [bl.source.shine()]
    outDict = {}
    for ii, (oe, screen) in enumerate(zip(bl.oes, bl.screens)):
        beam_s = screen.expose(beam=beam[ii])
        if 'slits' in oe.name:
            beam.append(oe.propagate(beam=beam[ii], needNewGlobal=True)[0])
        elif 'crl' in oe.name:
            beam.append(oe.multiple_refract(beam=beam[ii])[0])
        else:
            beam.append(oe.reflect(beam=beam[ii])[0])
        outDict[screen.name] = beam_s
            # note: two beam are output by the functions 'reflect'/'propagate'/...: 
            # [0] is the beam in global coordinates and [1] in local coordinates.
            # The global beam should be passed to subsequent element.
            # The local output can be suppressed by kwarg 'needLocal=False'
    screen = bl.screens[-1]
    beam_s = screen.expose(beam=beam[-1])
    outDict[screen.name] = beam_s
    return outDict
raycing.run.run_process = run_process


def define_plots(bl):
    plots = []
    for screen in bl.screens:
        saveName = ('./plots/' + subfolder + save_prefix + screen.name + '_y' 
                + str(int(screen.center[1])) + save_suffix + '.png')
        if 'slits' in screen.name:
            plot = XYCPlot(beam=screen.name,
                xaxis=xrtplot.XYCAxis(label="x",fwhmFormatStr='%.3f', limits=[-0.01,0.01]),
                yaxis=xrtplot.XYCAxis(label="z",fwhmFormatStr='%.3f', limits=[-0.05,0.05]),
                caxis=xrtplot.XYCAxis(label="energy",unit=r"eV",offset=E0,fwhmFormatStr='%.4f'),
                title=screen.name+ '_y' + str(int(screen.center[1])),
                saveName=str(saveName), aspect='auto')
        else:
            plot = XYCPlot(beam=screen.name,
                xaxis=xrtplot.XYCAxis(label="x",fwhmFormatStr='%.3f', limits=[-1.2,1.2]),
                yaxis=xrtplot.XYCAxis(label="z",fwhmFormatStr='%.3f'),
                caxis=xrtplot.XYCAxis(label="energy",unit=r"eV",offset=E0,fwhmFormatStr='%.4f'),
                title=screen.name+ '_y' + str(int(screen.center[1])),
                saveName=str(saveName))
            # plot = XYCPlot(beam=screen.name,
            #     xaxis=xrtplot.XYCAxis(label="x",fwhmFormatStr='%.3f'),
            #     yaxis=xrtplot.XYCAxis(label="z",fwhmFormatStr='%.3f'),
            #     caxis=xrtplot.XYCAxis(label="energy",unit=r"eV",offset=E0,fwhmFormatStr='%.4f'),
            #     title=screen.name+ '_y' + str(int(screen.center[1])),
            #     saveName=str(saveName))
        plots.append(plot)
        
        if do_plot=='all':
            saveName = ('./plots/' + subfolder + save_prefix + screen.name + '_y' 
                + str(int(screen.center[1])) + '_divVsE' + save_suffix + '.png')
            plot = XYCPlot(beam=screen.name,
                xaxis=xrtplot.XYCAxis(label="energy", unit="eV", offset=E0, limits=None),
                yaxis=xrtplot.XYCAxis(label="z'", unit="urad", factor=1E6, limits=None),
                caxis=xrtplot.XYCAxis(label="energy",unit="eV",offset=E0,fwhmFormatStr='%.4f',),
                title=screen.name+ '_y' + str(int(screen.center[1])),
                saveName=saveName, aspect='auto')
            plots.append(plot)

            saveName = ('./plots/' + subfolder + save_prefix + screen.name + '_y' 
                + str(int(screen.center[1])) + '_zVsE' + save_suffix + '.png')
            plot = XYCPlot(beam=screen.name,
                xaxis=xrtplot.XYCAxis(label="energy", unit="eV", offset=E0),
                yaxis=xrtplot.XYCAxis(label="z", unit="mm"),
                caxis=xrtplot.XYCAxis(label="energy",unit="eV",offset=E0,fwhmFormatStr='%.4f',),
                title=screen.name+ '_y' + str(int(screen.center[1])),
                saveName=saveName, aspect='auto')
            plots.append(plot)
    return plots


def main():
    bl = build_beamline()
    # E0 = bl.source.energies[0]
    bl.alignE = E0
    plots = define_plots(bl)
    xrtrun.run_ray_tracing(beamLine=bl, repeats=repeat, plots=plots)


if __name__ == '__main__':
    main()