import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
# plt.switch_backend('agg')

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')
import os

import re
from scipy.spatial.transform import Rotation

import xrt.backends.raycing.materials as rm


def get_perpendicular_vector(v):
    return np.cross([1,0,0], v)

def mono_kwargs(
    E0, 
    dE, 
    hkl, 
    miscut, 
    crystal_pitch_corr=[0,0,0,0],
    withLens=True, 
    screens=None
    ):
    
    Eratio = E0/18 # beam size scales as 1/E. Ref: 300um at 18keV
    FWHM_x = 1000
    FWHM_z = 300
    # FWHM_z = 300/Eratio

    screens_kwargs = []

    # SOURCE
    d_source = 200
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
        # 'polarization': 'h'
        'polarization': 'v'
        }
    angle = 0


    # CRYSTAL 0
    Si_kwargs = {
        'hkl': hkl,
        # 'tK': 300 # [K]
        'tK': 120 # [K]
        }
    Si_crystal = rm.CrystalSi(**Si_kwargs)
    Si_tth = 2*Si_crystal.get_Bragg_angle(E0)

    n = np.array([0, np.cos(angle), np.sin(angle)])
    center = [0,0,0]
    crystal0_kwargs = {
        'name': 'crystal0',
        'center': center,
        'pitch': 'auto',
        'positionRoll': 0,
        'material': Si_crystal,
        'alpha': 0
        }
    
    if (screens is None) or ('crystal0' in screens):
        screens_kwargs.append( {
            'name': 's_crystal0',
            'center': center,
            'z': get_perpendicular_vector(n)
            } )

    angle = angle + (Si_tth - 2*Si_crystal.get_dtheta(E0, 0) + crystal_pitch_corr[0])


    # CRYSTAL 1
    center = crystal0_kwargs['center'] + \
        np.array([0, np.cos(angle), np.sin(angle)])*200
    crystal1_kwargs = {
        'name': 'crystal1',
        'center': center,
        'pitch': 'auto',
        'positionRoll': np.pi,
        'material': Si_crystal,
        'alpha': miscut
        }

    if (screens is None) or ('crystal1' in screens):
        screens_kwargs.append( {
            'name': 's_crystal1',
            'center': center,
            'z': get_perpendicular_vector(n)
            } )

    angle = angle - (Si_tth - 2*Si_crystal.get_dtheta(E0, miscut) + crystal_pitch_corr[1])


    # CRL 1
    mBeryllium = rm.Material('Be', rho=1.848, kind='lens')
    q = 10000 # mm
    zmax = 0.2
    n = np.array([0, np.cos(angle), np.sin(angle)])
    center = crystal1_kwargs['center'] + n*q
    crl1_kwargs = {
        'name': 'crl1', 
        'center': center, 
        'material': mBeryllium,
        'zmax': zmax,
        'nCRL': (q, E0),
        }
    
    if (screens is None) or ('crl1' in screens):
        screens_kwargs.append( {
            'name': 's_crl1',
            'center': center,
            'z': get_perpendicular_vector(n)
            } )

    angle = angle+0


    # SLITS
    n = np.array([0, np.cos(angle), np.sin(angle)])
    center = crl1_kwargs['center'] + n*q
    slit_kwargs = {
        'name': 'slits',
        'center': center,
        'kind': ['left', 'right', 'bottom', 'top'],
        # 'opening': [-5, 5, -.0035, .0035]
        'opening': [-5, 5, -5, 5]
        }
    
    if (screens is None) or ('slits' in screens):
        screens_kwargs.append( {
            'name': 's_slits',
            'center': center,
            'z': get_perpendicular_vector(n)
            } )

    angle = angle+0


    # CRL 2
    n = np.array([0, np.cos(angle), np.sin(angle)])
    center = slit_kwargs['center'] + n*q
    crl2_kwargs = {
        'name': 'crl2', 
        'center': center, 
        'material': mBeryllium,
        'zmax': zmax,
        'nCRL': (q, E0),
        }
    
    if (screens is None) or ('crl2' in screens):
        screens_kwargs.append( {
            'name': 's_crl2',
            'center': center,
            'z': get_perpendicular_vector(n)
            } )

    angle = angle+0


    # CRYSTAL 2
    n = np.array([0, np.cos(angle), np.sin(angle)])
    center = crl2_kwargs['center'] + n*(q-200)
    crystal2_kwargs = {
        'name': 'crystal2',
        'center': center,
        'pitch': 'auto',
        'positionRoll': np.pi,
        'material': Si_crystal,
        'alpha': 0
        }

    if (screens is None) or ('crystal2' in screens):
        screens_kwargs.append( {
            'name': 's_crystal2',
            'center': center,
            'z': get_perpendicular_vector(n)
            } )

    angle = angle - (Si_tth - 2*Si_crystal.get_dtheta(E0, 0) + crystal_pitch_corr[2])


    # CRYSTAL 3
    n = np.array([0, np.cos(angle), np.sin(angle)])
    center = crystal2_kwargs['center'] + n*200
    crystal3_kwargs = {
        'name': 'crystal3',
        'center': center,
        'pitch': 'auto',
        'positionRoll': 0,
        'material': Si_crystal,
        'alpha': -miscut
        }
    
    if (screens is None) or ('crystal3' in screens):
        screens_kwargs.append( {
            'name': 's_crystal3',
            'center': center,
            'z': get_perpendicular_vector(n)
            } )

    angle = angle + (Si_tth - 2*Si_crystal.get_dtheta(E0, -miscut) + crystal_pitch_corr[3])


    # DETECTOR
    n = np.array([0, np.cos(angle), np.sin(angle)])
    center = crystal3_kwargs['center'] + n*300
    screens_kwargs.append( {
        'name': 'detector',
        'center': center,
        'z': get_perpendicular_vector(n)
        } )


    if withLens is True:
        oes_kwargs = [source_kwargs, crystal0_kwargs, crystal1_kwargs, crl1_kwargs, \
            slit_kwargs, crl2_kwargs, crystal2_kwargs, crystal3_kwargs]
    else:
        oes_kwargs = [source_kwargs, crystal0_kwargs, crystal1_kwargs, \
            slit_kwargs, crystal2_kwargs, crystal3_kwargs]
        screens_kwargs.pop(2)
        screens_kwargs.pop(3)
    
    return oes_kwargs, screens_kwargs