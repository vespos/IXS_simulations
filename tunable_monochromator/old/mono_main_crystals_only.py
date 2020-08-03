import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')

import re
from scipy.spatial.transform import Rotation

import xrt.backends.raycing as raycing
import xrt.plotter as xrtplot
import xrt.runner as xrtrun

from xrt.plotter import XYCAxis, XYCPlot
from xrt.backends.raycing.sources import GeometricSource
from xrt.backends.raycing.screens import Screen
from xrt.backends.raycing.oes import OE
from xrt.backends.raycing.materials import Crystal, CrystalSi

import fun


""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% """
""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% """
instrument = 'TXI'
E0 = 4400 # [eV]
dE = 0.0005 # [eV]
step = 0.05 # used for distE='lines' only
save_suffix = ''
n_refl = 4 # number of crystal to load (1-4)

# Source
source_kwargs = {
    'center': [0,0,0],
    'nrays': 1e5,
    # 'distE': 'normal', # lines, normal
    # 'energies': [E0,dE],
    'distE': 'lines',
    'energies': np.arange(E0-dE, E0+dE, step),
    'polarization': 'h'
    }


# Crystals, mirror and beam directions
miscut = np.radians(21)
d01 = 30 # distance between crystal 0 and 1 [mm]
d12 = 120 # distance between crystal 1 and 2 [mm]
Si_kwargs = {
    'hkl': [1,1,1],
    'tK': 300 # [K]
    }
Si_crystal = CrystalSi(**Si_kwargs)
Si_tth = 2*Si_crystal.get_Bragg_angle(E0)

# Prepare kwargs for all 4 OE crystal in one loop
d = [100, d01, d12, d01] # distance to previous crystal (source for the first one)
miscuts = [0, miscut, 0, miscut]
n= [np.array([0,1,0])] # unit vector indicating propagationd direction of the beam
beam_angle = [Si_tth, -Si_tth, -Si_tth, Si_tth] # angle that the beam make at each crystal
angle = np.cumsum(beam_angle) # beam angle in global frame
position_roll = [0, np.pi, np.pi, 0] # element facing up or down

xtal_pos = []
pitch = []
xtal_kwargs = []
for ii in range(n_refl):
    if ii==0:
        xtal_pos.append(n[0]*d[0])
    else:
        xtal_pos.append(xtal_pos[ii-1] + n[ii]*d[ii])
    n.append(np.array([0,
        np.cos(angle[ii]),
        np.sin(angle[ii])])) # unitary vector along the beam direction (global coord)
    pitch.append((Si_tth/2 - Si_crystal.get_dtheta(E0, alpha=miscuts[ii]) - miscuts[ii])
            *np.cos(position_roll[ii]))
            # the factor np.cos(position_roll) accounts for sign reversal of flipped elements
    xtal_kwargs.append({
        'name': 'SiCrystal'+str(ii),
        'center': xtal_pos[ii],
        # 'pitch': pitch[ii],
        'pitch': 'auto',
        'roll': 0,
        'yaw': 0,
        'positionRoll': position_roll[ii],
        'material': Si_crystal,
        'alpha': miscuts[ii]
        })
print(pitch)



# Screens
def get_perpendicular_vector(v):
    return np.cross([1,0,0], v)

s_kwargs = []
for ii, center in enumerate(xtal_pos):
    s_kwargs.append({
    'name': 's'+str(ii),
    'center': xtal_pos[ii],
    'z': get_perpendicular_vector(n[ii])
    })

detector_kwargs = {
    'name': 'detector',
    'center': s_kwargs[-1]['center']+[0,50,0],
    'z': get_perpendicular_vector(n[-1])
    }



""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% """
""" %%%%%%%%%% ELEMENTS CONSTRUCTORS %%%%%%%%%% """
""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% """

""" %%%%% SOURCE %%%%% """
class LCLS_source(GeometricSource):
    """
    LCLS source constructor.

    instance:
        Geometric source with LCLS instrument source parameters
    
    Note: by default the source is created at the sample position of the selected instrument.
    The instrument 'origin' (default) can be use to place the source at the real LCLS source 
    position
    """
    def __init__(self, bl=None, instrument='origin', d_offset=0, **kwargs):
        """
        inputs:
            bl: instance of the beamline class from xrt
            instrument: LCLS instrument ('TXI','XPP','XCS','MFX','CXI','MEC')
            d_offset: offset distance from the sample position
            **kwargs: kwargs argument to pass to GeometricSource instance
        """
        GeometricSource.__init__(self, bl=bl, name=instrument, **kwargs)
        FWHM, div = self.LCLS_source_properties()
        d = {'origin':0, 'TXI':140, 'XPP':150, 'XCS':400, 'MFX':415, 'CXI':425, 'MEC':440} #[m]
        # Sample position for each instrument
        d = d[instrument] + d_offset
        FWHM = FWHM + div*d #um
        self.dx = FWHM/1000/2.634 # sigma [mm]
        self.dz = FWHM/1000/2.634 # sigma [mm]
        self.dxprime = div /1e6 # [rad]
        self.dzprime = div /1e6 # [rad]

    def LCLS_source_properties(self):
        """
        LCLS source size and divergence (from Hasan's spreadsheet)
        input:
            en: photon energy [eV]
        output:
            FWHM: source size [mm]
            div: source divergence [urad]
        """
        E0 = self.energies[0]        
        if E0>=2000 and E0<8000:
            FWHM = 45*(5000/E0) #um
            div = 3.4*(5000/E0) #um
        if E0>=8000 and E0<=13000:
            FWHM = 37*(8000/E0)**0.25 #urad
            div = 2*(8000/E0)**0.25 #urad
        return FWHM, div




""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% """
""" %%%%%%%%%% XRT MAIN FUNCTIONS %%%%%%%%%% """
""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% """

def build_beamline():
    bl = raycing.BeamLine()

    # Source
    bl.source = LCLS_source(bl=bl, instrument=instrument, **source_kwargs)

    # Crystals
    xtals = []
    for kwargs in xtal_kwargs:
        xtals.append(OE(bl=bl, **kwargs))
    bl.xtals = xtals

    # Screens
    for kwargs in s_kwargs:
        Screen(bl=bl, **kwargs)
    bl.detector = Screen(bl=bl, **detector_kwargs)
    return bl


def run_process(bl):
    beam = [bl.source.shine()]
    outDict = {}
    for ii, xtal in enumerate(bl.xtals):
        beam_s = bl.screens[ii].expose(beam=beam[ii])
        beam.append(xtal.reflect(beam=beam[ii])[0])
        outDict[bl.screens[ii].name] = beam_s
            # note: two beam are output by the function 'reflect': 
            # [0] is the beam in global coordinates and [1] in local coordinates.
            # The global beam should be passed to subsequent element.
            # The local output can be suppressed by kwarg 'needLocal=False'
    beam_det = bl.detector.expose(beam=beam[-1])
    outDict['detector'] = beam_det
    return outDict
raycing.run.run_process = run_process


def define_plots(bl):
    plots = []
    for screen in bl.screens:
        saveName = ('./plots/' + screen.name + '_y' + str(int(screen.center[1]))
            + save_suffix + '.png')
        plot = XYCPlot(beam=screen.name,
            xaxis=xrtplot.XYCAxis(label="x"),
            yaxis=xrtplot.XYCAxis(label="z"),
            caxis=xrtplot.XYCAxis(label="energy",unit=r"eV", offset=E0),
            title=screen.name+ '_y' + str(int(screen.center[1])),
            saveName=str(saveName))
        plots.append(plot)
        
        saveName = ('./plots/' + screen.name + '_y' + str(int(screen.center[1]))
            + '_divVsE' + save_suffix + '.png')
        plot = XYCPlot(beam=screen.name,
            xaxis=xrtplot.XYCAxis(label="energy", unit="eV", offset=E0),
            yaxis=xrtplot.XYCAxis(label="z'", unit="rad"),
            caxis=xrtplot.XYCAxis(label="energy",unit="eV", offset=E0),
            aspect='auto',
            title=screen.name+ '_y' + str(int(screen.center[1])),
            saveName=saveName)
        plots.append(plot)
    return plots


def main():
    bl = build_beamline()
    # E0 = bl.source.energies[0]
    bl.alignE = E0
    plots = define_plots(bl)
    xrtrun.run_ray_tracing(beamLine=bl, repeats=1, plots=plots)


if __name__ == '__main__':
    main()