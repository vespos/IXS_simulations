import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')

import xrt.backends.raycing as raycing
import xrt.plotter as xrtplot
import xrt.runner as xrtrun

from xrt.plotter import XYCAxis, XYCPlot
from xrt.backends.raycing.sources import GeometricSource
from xrt.backends.raycing.screens import Screen
from xrt.backends.raycing.oes import OE, ParabolicalMirrorParam, EllipticalMirrorParam
from xrt.backends.raycing.materials import Material

def get_perpendicular_vector(v):
    return np.cross([1,0,0], v)

E0 = 4400
m_mat = Material(elements='Au', rho=19.32, kind='mirror')
m_pitch = np.radians(0.2) # mirror pitch
q = 500 # focal distance

D = 50E3 # distance from the source
angle = 2*m_pitch
n = np.array([0,
    np.cos(angle),
    np.sin(angle)]) # unit vector along the beam after first mirror


def build_beamline():
    bl = raycing.BeamLine()

    # Source
    bl.source = GeometricSource(
        bl=bl, dx=0.02, dz=0.02, dxprime=4E-6, dzprime=4E-6, energies=[E0])
    # bl.source = GeometricSource(
    #     bl=bl, dx=0.02, dz=0.02, dxprime=0, dzprime=0, energies=[E0])

    # Optical elements
    bl.oes = []
    ParabolicalMirrorParam(
        bl=bl,
        center=np.array([0,D,0]),
        pitch=m_pitch,
        q=q,
        material=m_mat,
        isCylindrical=True,
        limPhysY=[-200, 200]
        )

    ParabolicalMirrorParam(
        bl=bl,
        center=bl.oes[0].center+2*q*n,
        pitch=-m_pitch,
        positionRoll=np.pi,
        p=q,
        material=m_mat,
        isCylindrical=True,
        limPhysY=[-200, 200]
        )

    # Screens
    bl.s1 = Screen(
        bl=bl,
        name='source',
        center=np.asarray(bl.source.center)
        )
    bl.s2 = Screen(
        bl=bl,
        name='focus',
        center=bl.oes[0].center+q*n,
        z=get_perpendicular_vector(n)
        )
    bl.s3 = Screen(
        bl=bl,
        name='detector',
        center=bl.oes[1].center+np.array([0,300,0]),
        )

    print('\n\n')
    print('OEs position:')
    print(bl.oes[0].center)
    print(bl.oes[1].center)
    print('\n')
    print('\nScreens position:')
    print(bl.s1.center)
    print(bl.s2.center)
    print(bl.s3.center)
    print('\n\n')
    return bl


def run_process(bl):
    beam = [bl.source.shine()]
    for ii, oe in enumerate(bl.oes):
        beam.append(oe.reflect(beam=beam[ii])[0])
    outDict = {}
    outDict['source'] = bl.s1.expose(beam[0])
    outDict['focus'] = bl.s2.expose(beam=beam[1])
    outDict['detector'] = bl.s3.expose(beam=beam[-1])
    return outDict
raycing.run.run_process = run_process


def define_plots(bl):
    plots = []
    for screen in bl.screens:
        # x-z plot:
        saveName = './plots/dist_'+str(int(D/1000))+'_'+screen.name+'.png'
        plot = XYCPlot(beam=screen.name,
            xaxis=xrtplot.XYCAxis(label="x",fwhmFormatStr='%.3f'),
            yaxis=xrtplot.XYCAxis(label="z",fwhmFormatStr='%.3f'),
            caxis=xrtplot.XYCAxis(label="energy",unit=r"eV",offset=E0,fwhmFormatStr='%.4f'),
            title=screen.name+ '_y' + str(int(screen.center[1])),
            saveName=str(saveName))
        plots.append(plot)
        # E-z' plot:
        saveName = './plots/dist_'+str(int(D/1000))+'_'+screen.name+'._div.png'
        plot = XYCPlot(beam=screen.name,
                xaxis=xrtplot.XYCAxis(label="energy", unit="eV", offset=E0),
                yaxis=xrtplot.XYCAxis(label="z'", unit="urad", factor=1E6),
                caxis=xrtplot.XYCAxis(label="energy",unit="eV",offset=E0,fwhmFormatStr='%.4f',),
                title=screen.name+ '_y' + str(int(screen.center[1])),
                saveName=saveName, aspect='auto')
        plots.append(plot)
    return plots


def main():
    bl = build_beamline()
    bl.alignE = E0
    plots = define_plots(bl)
    xrtrun.run_ray_tracing(beamLine=bl, repeats=1, plots=plots)


if __name__ == '__main__':
    main()