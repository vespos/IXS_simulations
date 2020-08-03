import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sin, cos, tan, arctan, degrees, radians

sys.path.append(r'C:\Users\espov\Documents\Python\xrt')
sys.path.append('/mnt/c/Users/espov/Documents/Python/xrt')
import xrt.backends.raycing as raycing
import xrt.runner as xrtrun
import xrt.plotter as xrtplot
from xrt.backends.raycing.materials import Material
from xrt.backends.raycing.materials import Multilayer
from xrt.backends.raycing.screens import Screen

sys.path.append(r'C:\Users\espov\Documents\LCLS_projects\IXS\spectrometer')
sys.path.append('/mnt/c/Users/espov/Documents/LCLS_projects/IXS/spectrometer/')
from elements import LMultilayer
from elements import PMirror
from elements import LSource

import matplotlib as mpl
mpl.use('Agg')

import argparse
parser = argparse.ArgumentParser(description='Optimize pitch of crystals')
parser.add_argument('-a','--arg', type=float, help='argument')
parser.add_argument('-f','--file', type=str, help='file to save data')
args = parser.parse_args()
angle_i = args.arg
save_file = args.file
# ---------------------------------------------------------------------------

def get_perpendicular_vector(v):
    return np.cross([1,0,0], v)

def get_propagation_vector(angle):
    # unitary vector along the beam propagation direction
    return np.array([0, np.cos(angle), np.sin(angle)])

if not save_file:
    save_file = './test_reflectivity.dat'

if not angle_i:
    angle_i = 1.5 # grazing angle of incidence (deg)
beamInDotNormals = sin(radians(angle_i))

# source properties
nrays = 1E4
E0 = 8798
dE = 0.04
bsize = (1, 0.01)
# bdiv = (0, 0)
bdiv = (100E-6, 5E-3)
enstep = 0.005

# parabolic ML mirror
p = 0.35  # mm <- parabolic parameter, aka radius at vertex
S = 200  # mm <- distance between focal point and mirror center, focal distance
mirr_len = 150
mirr_width = 7
braggML = radians(angle_i)  # central 'bragg' angle of mirror #1
slope_error = 0  # in microradians

mSi = Material('Si', rho=2.33)
mRu = Material('Ru', rho=11.0)
mC = Material('C', rho=2.26)
mL = LMultilayer(tLayer=mC, tThickness=30, bLayer=mRu, bThickness=30,
                nPairs=150, substrate=mSi, p=p, E0=E0)

mirr_pos = [0, S, 0]
mirr_p = [0, 0, 0] #focal point is source
mirr_out = braggML*2

n0 = [0,1,0]
n1 = get_propagation_vector(mirr_out)

def build_beamline():
    bl = raycing.BeamLine()
    bl.source = LSource(bl, nrays=nrays, energies=[E0, dE], enstep=enstep,
                    distx='flat', dx=[-bsize[0]/2,bsize[0]/2],
                    distz='flat', dz=[-bsize[1]/2,bsize[1]/2],
                    distxprime='flat', dxprime=[-bdiv[0],bdiv[0]],
                    distzprime='flat', dzprime=[-bdiv[1],bdiv[1]])

    bl.mirr = PMirror(bl, center=mirr_pos, material=mL,
                    pitch=braggML,
                    limPhysY=[-mirr_len/2, mirr_len/2],
                    limPhysX=[-mirr_width/2, mirr_width/2],
                    f1=mirr_p, slope_error=slope_error)

    bl.s0 = Screen(bl, center=mirr_pos)
    bl.det = Screen(bl, center=mirr_pos+n1*100, z=get_perpendicular_vector(n1))
    return bl

def run_process(bl):
    beam_r0 = bl.source.shine()
    beam_r1 = bl.mirr.reflect(beam_r0)[0]

    outDict = {
        's0': bl.s0.expose(beam_r0),
        'det': bl.det.expose(beam_r1)
        }
    return outDict

def define_plots(bl):
    plots = []
    plot = xrtplot.XYCPlot(beam='s0',
        xaxis=xrtplot.XYCAxis(label="x",fwhmFormatStr='%.3f',limits=[-1,1]),
        yaxis=xrtplot.XYCAxis(label="z",fwhmFormatStr='%.3f', limits=[-1.5,1]),
        caxis=xrtplot.XYCAxis(label="energy",unit=r"eV",offset=E0,fwhmFormatStr='%.4f'))
    plots.append(plot)
    plot = xrtplot.XYCPlot(beam='det',
        xaxis=xrtplot.XYCAxis(label="x",fwhmFormatStr='%.3f',limits=[-1,1]),
        yaxis=xrtplot.XYCAxis(label="z",fwhmFormatStr='%.3f', limits=[-1.5,1]),
        caxis=xrtplot.XYCAxis(label="energy",unit=r"eV",offset=E0,fwhmFormatStr='%.4f'))
    plots.append(plot)
    return plots

raycing.run.run_process = run_process

if __name__ == '__main__':
    bl = build_beamline()
    bl.alignE = E0
    plots = define_plots(bl)
    xrtrun.run_ray_tracing(beamLine=bl, repeats=1, plots=plots)

    s_ind = -1 # screen from which the intensity should be taken
    intensity = plots[s_ind].intensity
    print('\n\nIntensity: ')
    print(intensity)
    with open(save_file, "a") as f:
        data = '{},{}\n'.format(str(angle_i), str(intensity))
        f.write(data)
    print('\n\n\n')