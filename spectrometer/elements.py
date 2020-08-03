
__author__ = "Hasan Yavas, Joel Bertinshaw, Hlynur Gretarsson"
__date__ = "2019/05"

import inspect
import numpy as np
from numpy import pi, sin, cos, tan, arcsin, arctan, degrees, radians, sqrt

import sys
sys.path.append(r"/Users/yavas/xrt-master")

import xrt.backends.raycing as raycing
from xrt.backends.raycing.oes import OE, ParabolicalMirrorParam
from xrt.backends.raycing.materials import Multilayer
from xrt.backends.raycing.sources import GeometricSource
from xrt.backends.raycing.sources_beams import Beam
from xrt.backends.raycing.physconsts import PI2


class LSource(GeometricSource):


    def __init__(self, *args, **kwargs):
        self.enstep = kwargs.pop('enstep', None)
        self.ensig = kwargs.pop('ensig', None)
        self.enmix = kwargs.pop('enmix', 1)
        GeometricSource.__init__(self, *args, **kwargs)
        self.distE = 'lines'


    def make_energy(self):

        def _in1d_tol(a,b,tol):
            d=np.abs(a-b[:,np.newaxis])
            return np.any(d<=tol, axis=0)

        nrays = self.nrays
        E0, dE = self.energies
        Esig, Estep, mix = self.ensig, self.enstep, self.enmix

        if not Estep:
            return np.random.normal(E0, dE, nrays)

        if isinstance(Estep,list):
            en = np.array(Estep) + E0
        else:
            en = np.arange(E0-dE, E0+dE+Estep, Estep)

        if not Esig:
            return np.array(en)[np.random.randint(len(en), size=nrays)]

        x = np.arange(E0-(dE+dE/2),E0+(dE+dE/2+0.00001),0.00001)
        y = np.zeros(x.shape)
        y[_in1d_tol(x,en,1E-5)] = 1
        xl = np.linspace(-100,100,len(x))

        func = (1-mix) * np.exp(-1*xl**2/(2*Esig**2))
        func += mix * Esig/(xl**2+Esig**2)

        y = np.convolve(y,func,mode='same')
        y /= y.sum()
        return np.random.choice(x,size=int(self.nrays),p=y)


    def make_polarization(self, bo):

        def _fill_beam(Jss, Jpp, Jsp, Es, Ep):
            bo.Jss.fill(Jss)
            bo.Jpp.fill(Jpp)
            bo.Jsp.fill(Jsp)
            if hasattr(bo, 'Es'):
                bo.Es.fill(Es)
                if isinstance(Ep, str):
                    bo.Ep[:] = np.random.uniform(size=int(nrays)) * 2**(-0.5)
                else:
                    bo.Ep.fill(Ep)

        polarization = self.polarization
        nrays = self.nrays
        if (polarization is None) or (polarization.startswith('un')):
            _fill_beam(0.5, 0.5, 0, 2**(-0.5), 'random phase')
        elif isinstance(polarization, tuple):
            if len(polarization) != 4:
                raise ValueError('wrong coherency matrix: must be a 4-tuple!')
            bo.Jss.fill(polarization[0])
            bo.Jpp.fill(polarization[1])
            bo.Jsp.fill(polarization[2] + 1j*polarization[3])
        else:
            if polarization.startswith('h'):
                _fill_beam(1, 0, 0, 1, 0)
            elif polarization.startswith('v'):
                _fill_beam(0, 1, 0, 0, 1)
            elif polarization == '+45':
                _fill_beam(0.5, 0.5, 0.5, 2**(-0.5), 2**(-0.5))
            elif polarization == '-45':
                _fill_beam(0.5, 0.5, -0.5, 2**(-0.5), -2**(-0.5))
            elif polarization.startswith('r'):
                _fill_beam(0.5, 0.5, 0.5j, 2**(-0.5), -1j * 2**(-0.5))
            elif polarization.startswith('l'):
                _fill_beam(0.5, 0.5, -0.5j, 2**(-0.5), 1j * 2**(-0.5))


    def shine(self, toGlobal=True, withAmplitudes=False, accuBeam=None):
        """.. Returned values: beamGlobal """

        try:
            self.bl._alignE = float(self.bl.alignE)
        except ValueError:
            self.bl._alignE = self.energies[0]

        bo = Beam(self.nrays)
        bo.state[:] = 1

        self.make_polarization(bo)
        self._apply_distribution(bo.y, self.disty, self.dy)
        self._apply_distribution(bo.x, self.distx, self.dx, bo)
        self._apply_distribution(bo.z, self.distz, self.dz, bo)
        self._apply_distribution(bo.a, self.distxprime, self.dxprime)
        self._apply_distribution(bo.c, self.distzprime, self.dzprime)

        ac = bo.a**2 + bo.c**2
        if sum(ac > 1) > 0:
            bo.b[:] = (ac + 1)**0.5
            bo.a[:] /= bo.b
            bo.c[:] /= bo.b
            bo.b[:] = 1.0 / bo.b
        else:
            bo.b[:] = (1 - ac)**0.5

        bo.E[:] = self.make_energy()

        if self.pitch or self.yaw:
            raycing.rotate_beam(bo, pitch=self.pitch, yaw=self.yaw)
        if toGlobal:
            raycing.virgin_local_to_global(self.bl, bo, self.center)
        raycing.append_to_flow(self.shine, [bo], inspect.currentframe())
        return bo


class LMultilayer(Multilayer):

    def __init__(self, *args, **kwargs):
        self.E0 = kwargs.pop('E0', None)
        self.p = kwargs.pop('p', None)
        self.wl = 12398.4191/self.E0
        Multilayer.__init__(self, *args, **kwargs)

        nt = self.tLayer.get_refractive_index(self.E0).real
        nb = self.bLayer.get_refractive_index(self.E0).real

        self.delta = 0.00012 #manually set to maximise flux, calculated = 0.000134

        self.delta = (nt-1)*self.tThicknessHigh + (nb-1)*self.bThicknessHigh  #correction to the refractive index is calculated
        self.delta = abs(self.delta) / self.d                                 #but it is not accurate, so close to the calculated value
                                                                                #should be tweaked, in this case 0.00012
    def dspacing(self, s):
        s = np.abs(s)
        th = pi/2 - arctan((2*s/self.p - 1)**0.5)
        d = self.wl / (2 * sin(th) * (1 - self.delta/(sin(th)**2)))
        return d

    def get_t_thickness(self, s, phi, iPair) :
        f = self.dspacing(s)
        return f/2

    def get_b_thickness(self, s, phi, iPair):
        f = self.dspacing(s)
        return f/2


class PMirror(ParabolicalMirrorParam):

    def __init__(self, *args, **kwargs):
        kwargs = self.__pop_kwargs(**kwargs)
        OE.__init__(self, *args, **kwargs)
        self.isParametric = True
        self.isCylindrical = True
        self.reset_pq(self.p, self.q, self.f1, self.f2)
        #maybe generate a profile
        #if self.slope_error:
        #    self.gen_slope_error()

    def __pop_kwargs(self, **kwargs):
        self.f1 = kwargs.pop('f1', None)
        self.f2 = kwargs.pop('f2', None)
        self.p = kwargs.pop('p', None)
        self.q = kwargs.pop('q', None)
        self.slope_error = kwargs.pop('slope_error',None)
        return kwargs

    def local_r_distorted(self, s, phi):
        """Distortion correction to the surface at (s, phi)
           (1) local_r += local_r_distorted
        """
        if not self.slope_error:
            return
        # height_error estimate using doi:10.1107/s1600577516007426
        # 1urad == 10nm, converted to mm
        height_error = self.slope_error * 10 * 1E-6
        return np.random.normal(scale=height_error, size=s.shape)

    def local_n_distorted(self, s, phi): #the reflected beam is shifted
        """Distortion to the local normal. Works in two possible ways:
           (1) Return d_pitch and d_roll rotation angles of the normal (i.e.
               rotations Rx and Ry). Works for parametric since the two
               rotations are around Cartesian axes and the local normal
               (local_n) is also a 3D vector in local xyz space.
           (2) local_n += local_n_distorted. Return a 3D vector that will be
               added to the local normal. The resulted vector will be
               normalized internally before calculating the reflected beam
               direction. A tuple of 3 arrays must be returned.
        """
        if not self.slope_error:
            return
        locx, locy, locz = self.param_to_xyz(s, phi, self.local_r(s, phi))
        slope_error = self.slope_error * 1E-6 #convert to radians
        r1 = 0 # slope error only along mirror
        #r1 = np.random.normal(scale=slope_error, size=locx.shape)
        r2 = np.random.normal(scale=slope_error, size=locy.shape)
        r3 = np.random.normal(scale=slope_error, size=locz.shape)
        return (r1,r2,r3)

    ## make a profile
    ## https://journals.iucr.org/s/issues/2016/04/00/ie5161/ie5161.pdf
    # def gen_slope_error(self):
    #     xmax, ymax = self.limPhysX[1], self.limPhysY[1]
    #     dx, dy = 0.1, 0.1
    #     x = np.arange(-xmax, xmax+dx, dx)
    #     y = np.arange(-ymax, ymax+dy, dy)
    #     #e.g. gaussian bump
    #     #z = 2.32e-4 * np.exp(-x[:, np.newaxis]**2/20**2 - y**2/150**2)
    #     a, b = np.gradient(z)
    #     a = np.arctan(a/dx)
    #     b = np.arctan(b/dy)
    #     self.warpNX, self.warpNY = len(x), len(y)
    #     self.warpZ = ndimage.spline_filter(z)
    #     self.warpA = ndimage.spline_filter(a)
    #     self.warpB = ndimage.spline_filter(b)
