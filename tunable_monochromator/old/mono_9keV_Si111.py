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
from xrt.backends.raycing.oes import OE, ParabolicalMirrorParam
from xrt.backends.raycing.materials import Material, Crystal, CrystalSi
from xrt.backends.raycing.apertures import RectangularAperture, DoubleSlit



""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% """
""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% """
instrument = 'TXI'
E0 = 9000 # [eV]
dE = 0.5 # [eV]
step = 0.12 # used for distE='lines' only
subfolder = '9keV/test_path/'
save_prefix = ''
save_suffix = ''
build_oe = [True, True, True, True, True, True, True] # decide whether an element is built (in order)
oes = ['crystal', 'crystal', 'parabolic_mirror', 'rectangular_aperture',
    'parabolic_mirror', 'crystal', 'crystal'] # oes category
# do_plot = ["z'", "z", "path"]
do_plot = ["z'"]

if not os.path.exists('./plots/'+subfolder):
    os.mkdir('./plots/'+subfolder)

def get_perpendicular_vector(v):
    return np.cross([1,0,0], v)

### Source ###
source_kwargs = {
    'center': [0,0,0],
    'd_offset': 2.5, # [m]
    'nrays': 5e6,
    # 'distE': 'normal', # lines, normal, flat
    # 'energies': [E0,dE],
    # 'distE': 'lines',
    # 'energies': np.arange(E0-dE, E0+dE, step),
    'distE': 'flat',
    'energies': [E0-dE, E0+dE],
    'polarization': 'h'
    }


### Optical elements ###

# crystals
miscut = np.radians(7.7)
# miscut = np.radians(10)
Si_kwargs = {
    'hkl': [1,1,1],
    'tK': 300 # [K]
    }
Si_crystal = CrystalSi(**Si_kwargs)
Si_tth = 2*Si_crystal.get_Bragg_angle(E0)

# mirrors
m_mat = Material(elements='Au', rho=19.32, kind='mirror')
m_pitch = np.radians(0.4)
f = 300 # focal distance

# slits
#  left, right, bottom, top]
# opening = [-1, 1, 0.009, .010]
opening = [-1, 1, -.5, .5]


# distances
ds0 = 100 # distance between source and crystal 0
d01 = 30 # distance between crystal 0 and 1 [mm]
d1m0 = f # distance between crystal 1 and mirror 0 [mm]
dm0slit = f # distance between mirror 0 and slit
dslitm1 = f # distance between slit and mirror 1 [mm]
dm12 = f-d01 # distance between mirror 1 and crystal 2 [mm]

# Prepare kwargs for all OE elements in one loop
d = [ds0, d01, d1m0, dm0slit, dslitm1, dm12, d01] # distance to previous element (absolute)
position_roll = [0, np.pi, 0, 0, np.pi, np.pi, 0] # element facing up or down
miscuts = [0, miscut, None, None, None, 0, -miscut] # None for non-crystal elements
# miscuts = [0, 0, None, None, 0, 0] # None for non-crystal elements
qp = [None, None, 'q', None, 'p', None, None] # None for non-curved mirror elements

beam_angle = np.array(
    [Si_tth, -Si_tth, 2*m_pitch, 0, -2*m_pitch, -Si_tth, Si_tth]) # angle that the beam make at each oe
ttheta_corr = []
for misc in miscuts:
    if misc is not None:
        ttheta_corr.append( Si_crystal.get_dtheta(E0, misc) )
    else:
        ttheta_corr.append(0)
# manual_corr = [0.1E-6, 371.6E-6, 0, 0, 254E-6, 232.8E-6]
# beam_angle = beam_angle - ttheta_corr #- manual_corr
angle = np.cumsum(beam_angle) # beam angle in (y-z) global frame afer each element (0 = parallel to y)

n= [np.array([0.,1.,0.])] # unit vector indicating propagation of the beam. Starts along y (source)
oes_pos = []
pitch = []
oes_kwargs = []
s_kwargs = []
for ii, oe in enumerate(oes):
    if build_oe[ii] is False:
        beam_angle[ii] = 0
        if ii==0:
            oes_pos.append(n[0]*d[0])
        else:
            oes_pos.append(oes_pos[ii-1] + n[ii]*d[ii])
        n.append(n[ii]) # unitary vector along the beam direction (global coord)
        pitch.append(0)
        continue

    if ii==0:
        oes_pos.append(n[0]*d[0])
    else:
        oes_pos.append(oes_pos[ii-1] + n[ii]*d[ii])
    
    n.append(np.array([0,
        np.cos(angle[ii]),
        np.sin(angle[ii])])) # unitary vector along the beam direction (global coord)
    
    if oe=='crystal':
        pitch.append((Si_tth/2 - Si_crystal.get_dtheta(E0, alpha=miscuts[ii]) - miscuts[ii])
            *np.cos(position_roll[ii]))
            # the factor np.cos(position_roll) accounts for sign reversal of flipped elements
        oes_kwargs.append({
            'name': oe+str(ii),
            'center': oes_pos[ii],
            # 'pitch': pitch[ii],
            'pitch': 'auto',
            'roll': 0,
            'yaw': 0,
            'positionRoll': position_roll[ii],
            'material': Si_crystal,
            'alpha': miscuts[ii]
            })
    
    elif oe=='parabolic_mirror':
        pitch.append(m_pitch*np.cos(position_roll[ii]))
            # the factor np.cos(position_roll) accounts for sign reversal of flipped elements
        oes_kwargs.append({
            'name': oe+str(ii),
            'center': oes_pos[ii],
            qp[ii]: f,
            'pitch': pitch[ii],
            'roll': 0,
            'yaw': 0,
            'positionRoll': position_roll[ii],
            'material': m_mat,
            'isCylindrical': True,
            'limPhysY': [-50, 50]
            })

    elif oe=='flat_mirror':
        pitch.append(m_pitch*np.cos(position_roll[ii]))
            # the factor np.cos(position_roll) accounts for sign reversal of flipped elements
        oes_kwargs.append({
            'name': oe+str(ii),
            'center': oes_pos[ii],
            'pitch': pitch[ii],
            'roll': 0,
            'yaw': 0,
            'positionRoll': position_roll[ii],
            'material': m_mat
            })
    
    elif oe=='rectangular_aperture':
        r = Rotation.from_rotvec(np.array([1,0,0])*(m_pitch-np.pi/2)) # rotation to tilt slit
        pitch.append(0)
        oes_kwargs.append({
            'name': oe+str(ii),
            'center': oes_pos[ii],
            'kind': ['left', 'right', 'bottom', 'top'],
            'opening': opening,
            # 'opening': [-.5, .5, .041, .042],
            # 'z': r.apply(get_perpendicular_vector(n[3]))
            })

    else:
        print('Optical element category not implemented')
    
    # create a screen at each oe position that images the beam from the previous oe
    s_kwargs.append({
        'name': 's'+str(ii),
        'center': oes_pos[ii],
        'z': get_perpendicular_vector(n[ii])
        })

print('pitch:')
print(pitch)
print('\nbeam direction:')
print('\t'+'\n\t'.join([str(ii) for ii in n])+'\n')


# Special screens added manually
r = Rotation.from_rotvec(np.array([1,0,0])*(m_pitch-np.pi/2)) # rotation to tilt screen
focus_kwargs = {
    'name': 'on_focus',
    'center': s_kwargs[2]['center']+n[3]*f,
    'z': get_perpendicular_vector(n[3])
    }

focus_tilt_kwargs = {
    'name': 'on_focus_tilt',
    'center': s_kwargs[2]['center']+n[3]*f,
    'z': r.apply(get_perpendicular_vector(n[3]))
    }

detector_kwargs = {
    'name': 'detector',
    'center': s_kwargs[-1]['center']+n[-1]*300,
    'z': get_perpendicular_vector(n[-1])
    }

special_screens_kwargs = [focus_kwargs, focus_tilt_kwargs, detector_kwargs]
# special_screens_kwargs = [detector_kwargs]



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
        if 'd_offset' in kwargs:
            d_offset = kwargs.pop('d_offset')
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
        # div = 0
        return FWHM, div




""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% """
""" %%%%%%%%%% XRT MAIN FUNCTIONS %%%%%%%%%% """
""" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% """

def build_beamline():
    bl = raycing.BeamLine()

    # Source
    # bl.source = LCLS_source(bl=bl, instrument=instrument, **source_kwargs)
    pop = source_kwargs.pop('d_offset')
    # bl.source = GeometricSource(bl=bl, dx=100e-6, dz=100e-6, dxprime=2e-6, \
    #     dzprime=2e-6, **source_kwargs)
    bl.source = GeometricSource(bl=bl, dx=100e-6, dz=100e-6, dxprime=0, \
        dzprime=0, **source_kwargs)

    # Optical elements
    oes = []
    print('OEs position:')
    for kwargs in oes_kwargs:
        if 'crystal' in kwargs['name']:
            oes.append(OE(bl=bl, **kwargs))
        elif 'parabolic_mirror' in kwargs['name']:
            oes.append(ParabolicalMirrorParam(bl=bl, **kwargs))
        elif 'flat_mirror' in kwargs['name']:
            oes.append(OE(bl=bl, **kwargs))
        elif 'rectangular_aperture' in kwargs['name']:
            oes.append(RectangularAperture(bl=bl, **kwargs))
        else:
            print('OE not implemented (build beamline)')
        print('\t' + oes[-1].name + ': ' + str(oes[-1].center))
    bl.oes = oes

    # Screens
    print('\nScreens position:')
    for kwargs in s_kwargs:
        Screen(bl=bl, **kwargs)
        print('\t' + str(bl.screens[-1].center))
    for kwargs in special_screens_kwargs:
        Screen(bl=bl, **kwargs)
        print('\t' + str(bl.screens[-1].name) + ': ' + str(bl.screens[-1].center))

    # bl.on_focus = Screen(bl=bl, **focus_kwargs)
    # print('Focus:' + str(bl.screens[-1].center))
    # bl.detector = Screen(bl=bl, **detector_kwargs)
    # print('Detector:' + str(bl.screens[-1].center))
    print('\n')
    return bl


def run_process(bl):
    beam = [bl.source.shine()]
    outDict = {}
    for ii, (oe, screen) in enumerate(zip(bl.oes, bl.screens)):
        beam_s = screen.expose(beam=beam[ii])
        if 'aperture' in oe.name:
            beam.append(oe.propagate(beam=beam[ii], needNewGlobal=True)[0])
            # beam.append(oe.propagate(beam=beam[ii]))
        else:
            beam.append(oe.reflect(beam=beam[ii])[0])
        outDict[screen.name] = beam_s
            # note: two beam are output by the function 'reflect': 
            # [0] is the beam in global coordinates and [1] in local coordinates.
            # The global beam should be passed to subsequent element.
            # The local output can be suppressed by kwarg 'needLocal=False'
    for ii, s_kwargs in enumerate(special_screens_kwargs):
        screen = next((screen for screen in bl.screens if screen.name == s_kwargs['name']), None)
            # replace names that do not correspond to s_kwargs['name'] by None
        if 'focus' in screen.name:
            beam_s = screen.expose(beam=beam[3])
        elif 'detector' in screen.name:
            beam_s = screen.expose(beam=beam[-1])
        outDict[screen.name] = beam_s
    return outDict
raycing.run.run_process = run_process


def define_plots(bl):
    plots = []
    for screen in bl.screens:
        saveName = ('./plots/' + subfolder + save_prefix + screen.name + '_y' 
                + str(int(screen.center[1])) + save_suffix + '.png')
        if 'focus' in screen.name:
            plot = XYCPlot(beam=screen.name,
                xaxis=xrtplot.XYCAxis(label="x",fwhmFormatStr='%.3f'),
                yaxis=xrtplot.XYCAxis(label="z",fwhmFormatStr='%.3f'),
                caxis=xrtplot.XYCAxis(label="energy",unit=r"eV",offset=E0,fwhmFormatStr='%.4f'),
                title=screen.name+ '_y' + str(int(screen.center[1])),
                saveName=str(saveName), aspect='auto')
        else:
            plot = XYCPlot(beam=screen.name,
                xaxis=xrtplot.XYCAxis(label="x",fwhmFormatStr='%.3f',limits=[-1,1]),
                yaxis=xrtplot.XYCAxis(label="z",fwhmFormatStr='%.3f', limits=[-1.5,1]),
                caxis=xrtplot.XYCAxis(label="energy",unit=r"eV",offset=E0,fwhmFormatStr='%.4f'),
                title=screen.name+ '_y' + str(int(screen.center[1])),
                saveName=str(saveName))
        plots.append(plot)
        
        if 'z\'' in do_plot:
            saveName = ('./plots/' + subfolder + save_prefix + screen.name + '_y' 
                + str(int(screen.center[1])) + '_divVsE' + save_suffix + '.png')
            plot = XYCPlot(beam=screen.name,
                xaxis=xrtplot.XYCAxis(label="energy", unit="eV", offset=E0, limits=None),
                yaxis=xrtplot.XYCAxis(label="z'", unit="urad", factor=1E6, limits=None),
                caxis=xrtplot.XYCAxis(label="energy",unit="eV",offset=E0,fwhmFormatStr='%.4f',),
                title=screen.name+ '_y' + str(int(screen.center[1])),
                saveName=saveName, aspect='auto')
            plots.append(plot)

        if 'z' in do_plot:
            saveName = ('./plots/' + subfolder + save_prefix + screen.name + '_y' 
                + str(int(screen.center[1])) + '_zVsE' + save_suffix + '.png')
            plot = XYCPlot(beam=screen.name,
                xaxis=xrtplot.XYCAxis(label="energy", unit="eV", offset=E0),
                yaxis=xrtplot.XYCAxis(label="z", unit="mm"),
                caxis=xrtplot.XYCAxis(label="energy",unit="eV",offset=E0,fwhmFormatStr='%.4f',),
                title=screen.name+ '_y' + str(int(screen.center[1])),
                saveName=saveName, aspect='auto')
            plots.append(plot)

        if 'path' in do_plot:
            saveName = ('./plots/' + subfolder + save_prefix + screen.name + '_y' 
                + str(int(screen.center[1])) + '_path' + save_suffix + '.png')
            plot = XYCPlot(beam=screen.name,
                xaxis=xrtplot.XYCAxis(label="energy", unit="eV", offset=E0),
                yaxis=xrtplot.XYCAxis(label="path", unit="mm"),
                caxis=xrtplot.XYCAxis(label="energy",unit="eV",offset=E0,fwhmFormatStr='%.4f',),
                title=screen.name+ '_y' + str(int(screen.center[1])),
                saveName=saveName, aspect='auto')
    return plots


def main():
    bl = build_beamline()
    # E0 = bl.source.energies[0]
    bl.alignE = E0
    plots = define_plots(bl)
    xrtrun.run_ray_tracing(beamLine=bl, repeats=1, plots=plots)


if __name__ == '__main__':
    main()