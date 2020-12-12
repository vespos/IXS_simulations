import numpy as np

import sys
sys.path.append(r'C:\Users\espov\Documents\Python\xrt')

from xrt.backends.raycing.sources import GeometricSource

class LCLS_source(GeometricSource):
    """
    LCLS source constructor.

    instance:
        Geometric source with LCLS instrument source parameters
    
    Note: by default the source is created at the sample position of the selected instrument.
    The instrument 'origin' (default) can be use to place the source at the real LCLS source 
    position
    """
    def __init__(self, 
        bl=None, 
        instrument='origin', 
        E0=12000,
        d_offset=0, 
        **kwargs):
        """
        inputs:
            bl: instance of the beamline class from xrt
            instrument: LCLS instrument ('TXI','XPP','XCS','MFX','CXI','MEC')
            d_offset: offset distance from the sample position
            **kwargs: kwargs argument to pass to GeometricSource instance
        """
        GeometricSource.__init__(self, bl=bl, name=instrument, 
                energies=(E0,), **kwargs)
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