
__author__ = "Hasan Yavas, Joel Bertinshaw, Hlynur Gretarsson"
__date__ = "2019/05"

## only save figures
#import matplotlib as mpl
#mpl.use('Agg')

import sys
import pickle

import numpy as np
import matplotlib.pyplot as plt

import argparse

sys.path.append(r'C:\Users\espov\Documents\Python\xrt')

import xrt.backends.raycing as raycing
import xrt.plotter as xrtplot
import xrt.runner as xrtrun

from numpy import pi, sin, cos, tan, arctan, degrees, radians
from xrt.backends.raycing.oes import OE, ParabolicalMirrorParam
from xrt.backends.raycing.apertures import RectangularAperture
from xrt.backends.raycing.materials import Material, CrystalFromCell, CrystalSi
from xrt.backends.raycing.screens import Screen

from elements import LMultilayer, PMirror, LSource
import materials as materials


# details
fig_dir = './fig/'
dat_dir = './dat/'
title = '9keVv1'
nrays = 1E5
repeats = 10

#plot spectrometer in 3D
glow3D = False

#just show angle calculations
just_calc = False

# source
# E0 = 9000.0
E0 = 8798
# dE = 0.04
dE = 0.04
# bsize = (0.1, 0.01)  # 100 micron (H) by 30 micron (V) beam
bsize = (1, 0.01)
# bdiv = (5E-3, 5E-3)  # divergence is equal to that of 100mm analyzer at 10m
bdiv = (100E-6, 5E-3)

#enstep = 0 # gaussian profile same as GeometricSource
enstep = 0.006  # stepwidth for lines within dE range
#enstep = [0] # single value
#enstep = [0,0.005] # specific values

ensig = 0 # sigma for lines, 1.277 --> 3meV
enmix = 0 # basic pv mixing 0-gaussian, 1-lorentzian


# Detector / figures

# ss = ['S0','S1','S2','S3','S4']
#mode = ['X','E','Zp','EZp','D'] #XvsZ, EvsZ, Z'vsZ, EvsZ', detector

ss = ['S0','S1','S2','S3','S4']  # screen number after the number of reflections
# ss = ['S2', 'S3']
mode = ['D']
fig_bin = 256 #for standard figures

#Andor iKon L
det_pix = 512
det_pitch = 50 * 1E-3
det_dim = det_pitch * det_pix/2


# parabolic mirror, from specifications
p = 0.35  # mm <- parabolic parameter, aka radius at vertex
S = 200  # mm <- distance between focal point and mirror center, focal distance
mirr_len = 150
mirr_width = 7
braggML = radians(1.694)  # central 'bragg' angle of mirror #1
slope_error = 2  # in microradians

FF = 400 # mm <- second mirror to detector
# FF = 400*1.5 # mm <- second mirror to detector


mSi = Material('Si', rho=2.33)
mNi = Material('Ni', rho=8.9)
mRu = Material('Ru', rho=11.0)
mC = Material('C', rho=2.26)
mL = LMultilayer(tLayer=mC, tThickness=30, bLayer=mRu, bThickness=30,
                 nPairs=150, substrate=mSi, p=p, E0=E0)

# xtals
#crystalSi = CrystalSi(geom="Bragg", hkl=(4, 0, 0))

# crystalGe = materials.Ge400_300K
crystalGe = materials.Ge800_300K
#crystalGe = CrystalFromCell('Germanium', (4, 0, 0), a=5.657, atoms=['Ge']*8,
#                            atomsXYZ=[[0.0, 0.0, 0.0],[0.0, 0.5, 0.5],
#                                      [0.5, 0.0, 0.5],[0.5, 0.5, 0.0],
#                                      [.25, .25, .25],[.25, .75, .75],
#                                      [.75, .25, .75],[.75, .75, .25]])

alphaGe = radians(-80.0)
# alphaGe = radians(-81.0)
braggGe = crystalGe.get_Bragg_angle(E0) - crystalGe.get_dtheta(E0, alphaGe)

# braggGe = crystalGe.get_Bragg_angle(E0) - crystalGe.get_dtheta(E0, 0)
b1 = -sin(braggGe+alphaGe) / sin(braggGe-alphaGe)
b2 = 1/b1
#alphaGe = arctan(tan(braggGe)*(-b2-1)/(-b2+1))
alphaGe2=arctan(tan(braggGe)*(-b2-1)/(-b2+1))



# element angles
mirr1_p = [0, 0, 0] #focal point is source
mirr1_out = braggML*2

xtal1_pitch = mirr1_out + braggGe + alphaGe
xtal1_out = mirr1_out + braggGe*2 + crystalGe.get_dtheta(E0, alphaGe)
# xtal1_out += 441.7e-6 #manually optimised

xtal2_pitch = xtal1_out + braggGe + alphaGe2
xtal2_out = xtal1_out + braggGe*2  + crystalGe.get_dtheta(E0, alphaGe2)
xtal2_pitch -= radians(0.0047) # manually optimised if energy is off
# xtal2_pitch -= radians(0.007) # manually optimised if energy is off
# xtal2_out -= 531e-6 # manually optimised

braggML2 = radians(1.1985)  # central 'bragg' angle of mirror #2
mirr2_out = braggML2*2 + xtal2_out


# element spacing
xtal1_h = 200
xtal2_v = 25
mirr2_h = 200

#position of every elements should be calculated manually

mirr1_pos = [0, S, 0]
xtal1_pos = [0, xtal1_h + mirr1_pos[1], xtal1_h * tan(mirr1_out)]
xtal2_pos = [0, xtal2_v/tan(xtal1_out) + xtal1_pos[1], xtal2_v + xtal1_pos[2]]
mirr2_pos = [0, xtal2_pos[1]+mirr2_h, xtal2_pos[2] + mirr2_h*tan(xtal2_out)]

# focal point of montel determined by mirr2 position
Sfc = FF  # focal distance of the focusing mirror
mirr2_q = [0, mirr2_pos[1] - Sfc*cos(pi-mirr2_out),
              mirr2_pos[2] + Sfc*sin(pi-mirr2_out)]

det_mirr2 = Sfc #detector at focal point, distance can be tweaked
det_ang_offset = radians(-88.5) #detector angle offset from flat, in radians

slit_pos = [0, mirr2_pos[1] - (det_mirr2-20)*cos(pi-mirr2_out),
            mirr2_pos[2] + (det_mirr2-20)*sin(pi-mirr2_out)]


""" propagation direction """
angles = [0, mirr1_out, xtal1_out, xtal2_out, mirr2_out]
n = []
for angle in angles:
    n.append(np.array([0,
        np.cos(angle),
        np.sin(angle)])) # unitary vector along the beam direction (global coord)


## Setup simulation

def build_beamline():

    ## beamline
    BL = raycing.BeamLine()

    BL.source = LSource(BL, nrays=nrays, energies=[E0, dE],
                        enstep=enstep, ensig=ensig, enmix=enmix,
                        distx='flat', dx=[-bsize[0]/2,bsize[0]/2],
                        distz='flat', dz=[-bsize[1]/2,bsize[1]/2],
                        distxprime='flat', dxprime=[-bdiv[0],bdiv[0]],
                        distzprime='flat', dzprime=[-bdiv[1],bdiv[1]])

    BL.mirr1 = PMirror(BL, center=mirr1_pos, material=mL,
                       pitch=braggML,
                       limPhysY=[-mirr_len/2, mirr_len/2],
                       limPhysX=[-mirr_width/2, mirr_width/2],
                       f1=mirr1_p, slope_error=slope_error)

    BL.xtal1 = OE(BL, center=xtal1_pos, material=crystalGe,
                  pitch=xtal1_pitch, alpha=alphaGe)

    BL.xtal2 = OE(BL, center=xtal2_pos, material=crystalGe,
                  pitch=xtal2_pitch, alpha=alphaGe2)

    BL.mirr2 = PMirror(BL, center=mirr2_pos, material=mL,
                       pitch=braggML2, extraPitch=xtal2_out,
                       limPhysY=[-mirr_len/2, mirr_len/2],
                       limPhysX=[-mirr_width/2, mirr_width/2],
                       f2=mirr2_q, slope_error=slope_error)

    BL.slt1 = RectangularAperture(BL, center=slit_pos,
                                  opening=[-10.2, 10.2, -10, 10])

    ## virtual screens
    BL.s0 = Screen(BL, center=mirr1_pos)
    BL.s1 = Screen(BL, center=xtal1_pos, z=[0, -sin(mirr1_out), cos(mirr1_out)])
    BL.s2 = Screen(BL, center=xtal2_pos, z=[0, -sin(xtal1_out), cos(xtal1_out)])
    BL.s3 = Screen(BL, center=mirr2_pos, z=[0, -sin(xtal2_out), cos(xtal2_out)])
    BL.sdebug = Screen(BL, center=xtal2_pos+n[3]*15, z=[0, -sin(xtal2_out), cos(xtal2_out)])

    ## detector screen, generally placed at 2nd montel focal point
    s4_pos = [0, mirr2_pos[1] - det_mirr2*cos(pi-mirr2_out),
                 mirr2_pos[2] + det_mirr2*sin(pi-mirr2_out)]

    s4_ang = [0, -sin(mirr2_out + det_ang_offset),
                  cos(mirr2_out + det_ang_offset)]

    BL.s4 = Screen(BL, center=s4_pos, z=s4_ang)

    return BL


def run_process(BL):

    beam_r0 = BL.source.shine()
    beam_r1 = BL.mirr1.reflect(beam_r0)[0]
    beam_r2 = BL.xtal1.reflect(beam_r1)[0]
    beam_r3 = BL.xtal2.reflect(beam_r2)[0]
    beam_r4 = BL.mirr2.reflect(beam_r3)[0]
    beam_r5 = BL.slt1.propagate(beam_r4)



    outDict = {
               'S0': BL.s0.expose(beam_r0),
               'S1': BL.s1.expose(beam_r1),
               'S2': BL.s2.expose(beam_r2),
               'S3': BL.s3.expose(beam_r3),
               'S4': BL.s4.expose(beam_r4),
               'Sdebug': BL.sdebug.expose(beam_r3)
              }

    if glow3D:
        BL.prepare_flow()

    return outDict


def beamline_stats(BL, andor_a=det_ang_offset, andor_d=det_mirr2):

    mirr1_y, mirr1_z = BL.mirr1.center[1], BL.mirr1.center[2]
    xtal1_y, xtal1_z = BL.xtal1.center[1], BL.xtal1.center[2]
    xtal2_y, xtal2_z = BL.xtal2.center[1], BL.xtal2.center[2]
    mirr2_y, mirr2_z = BL.mirr2.center[1], BL.mirr2.center[2]
    slit1_y, slit1_z = BL.slt1.center[1],  BL.slt1.center[2]
    andor_y, andor_z = BL.s4.center[1], BL.s4.center[2]
    mirr1_a, mirr2_a = BL.mirr1.pitch, BL.mirr2.pitch+BL.mirr2.extraPitch
    xta11_a, xtal2_a = BL.xtal1.pitch, BL.xtal2.pitch

    disp = [f'mirr1_y:    {mirr1_y:12.6f}   mirr1_z:    {mirr1_z:12.6f}',
            f'xtal1_y:    {xtal1_y:12.6f}   xtal1_z:    {xtal1_z:12.6f}',
            f'xtal2_y:    {xtal2_y:12.6f}   xtal2_z:    {xtal2_z:12.6f}',
            f'mirr2_y:    {mirr2_y:12.6f}   mirr2_z:    {mirr2_z:12.6f}',
             '',
            f'slit1_y:    {slit1_y:12.6f}   slit1_z:    {slit1_z:12.6f}',
            f'andor_y:    {andor_y:12.6f}   andor_z:    {andor_z:12.6f}',
            f'andor_dist: {andor_d:12.6f}   andor_pitch:{andor_a:12.8f}',
             '',
            f'mirr1_pitch:{mirr1_a:12.8f}   mirr1_out:  {mirr1_out:12.8f}',
            f'xtal1_pitch:{xta11_a:12.8f}   xtal1_out:  {xtal1_out:12.8f}',
            f'xtal2_pitch:{xtal2_a:12.8f}   xtal2_out:  {xtal2_out:12.8f}',
            f'mirr2_pitch:{mirr2_a:12.8f}   mirr2_out:  {mirr2_out:12.8f}',
             '',
            f'xtal1_bragg:{braggGe:12.8f}   xtal2_bragg:{braggGe:12.8f}',
            f'xtal1_alpha:{alphaGe:12.8f}   xtal2_alpha:{alphaGe2:12.8f}',
            f'xtal1_b:    {b1:12.8f}   xtal2_b:    {b2:10.8f}']
    disp = '\n'.join(disp)

#    stats = {'mirr1_y': mirr1_y, 'mirr1_z': mirr1_z, 'mirr1_pitch': mirr1_a,
#              'xtal1_y': xtal1_y, 'xtal1_z': xtal1_z ,'xta11_pitch': xta11_a,
#              'xtal2_y': xtal2_y, 'xtal2_z': xtal2_z, 'xtal2_pitch': xtal2_a,
#              'mirr2_y': mirr2_y, 'mirr2_z': mirr2_z, 'mirr2_pitch': mirr2_a,
#              'andor_y': andor_y, 'andor_z': andor_z, 'andor_pitch': andor_a,
#              'mirr1_out': mirr1_out, 'xtal1_out': xtal1_out,
#              'xtal2_out': xtal2_out, 'mirr2_out': mirr2_out, 'andor_d': andor_d,
#              'xtal1_bragg':braggGe, 'xtal1_alpha':alphaGe, 'xtal1_b':b1,
#              'xtal2_bragg':braggGe, 'xtal2_alpha':alphaGe, 'xtal2_b':b2,
#              'E0':E0, 'dE':dE, 'bsize':bsize, 'bdiv':bdiv, 'source': source,
#              'Estep': Estep, 'nrays': nrays, 'repeats': repeats, 'disp':disp}

    print(disp+'\n')
    return disp


def define_plots(screen_names, mode=['X','E','Zp','EZp','D'],
                 step=None, step2=None):

    plots = []
    if step is not None and step2 is not None:
        save_ext = f'_{step:.3f}_{step2:.3f}.png'
    elif step is not None:
        save_ext = f'_{step:.3f}.png'
    else:
        save_ext = '.png'

    for od in screen_names:
        if 'X' in mode:
            plots.append(xrtplot.XYCPlot(od,aspect='auto',
                         yaxis=xrtplot.XYCAxis("z", "mm", bins=fig_bin),
                         xaxis=xrtplot.XYCAxis("x", "mm", bins=fig_bin),
                         saveName=f'{fig_dir}{title}_X_{od}{save_ext}'))
        if 'E' in mode:
            plots.append(xrtplot.XYCPlot(od,aspect='auto',
                         yaxis=xrtplot.XYCAxis("z", "mm", bins=fig_bin),
                         xaxis=xrtplot.XYCAxis("energy", "eV", bins=fig_bin),
                         saveName=f'{fig_dir}{title}_E_{od}{save_ext}'))

        if 'Zp' in mode:
            plots.append(xrtplot.XYCPlot(od,aspect='auto',
                         yaxis=xrtplot.XYCAxis("z", "mm", bins=fig_bin),
                         xaxis=xrtplot.XYCAxis("z'", "µrad", bins=fig_bin),
                         saveName=f'{fig_dir}{title}_Zp_{od}{save_ext}'))

        if 'EZp' in mode:
            plots.append(xrtplot.XYCPlot(od,aspect='auto',
                         yaxis=xrtplot.XYCAxis("energy", "eV", bins=fig_bin),
                         xaxis=xrtplot.XYCAxis("z'", "µrad", bins=fig_bin),
                         saveName=f'{fig_dir}{title}_EZp_{od}{save_ext}'))

        if 'D' in mode:
            plots.append(xrtplot.XYCPlot(od,aspect='auto',
                         yaxis=xrtplot.XYCAxis("z", "mm", ppb=1, bins=det_pix,
#                                               limits=[-det_dim, det_dim]),
                                               limits=[-5.15, 6]),
                         xaxis=xrtplot.XYCAxis("x", "mm", ppb=1, bins=det_pix,
#                                               limits=[-det_dim, det_dim]),
                                               limits=[-7, 7]),
                         saveName=f'{fig_dir}{title}_D_{od}{save_ext}'))

    for plot in plots:

        sn = plot.saveName[:-4].split('/')[-1]
        od = sn.split('_')[2]
        pm = sn.split('_')[1]

        if od == 'S0':
            if pm == 'E':
                plot.xaxis.limits = [E0-0.01,E0+0.01]
                plot.xaxis.offset = E0
        if od == 'S2':
            if pm == 'EZp':
                plot.yaxis.limits = [E0-dE,E0+dE]
                plot.xaxis.limits = [-800, 800]
            else:
                plot.yaxis.limits = [-60, 60]
            if pm == 'X':
               plot.xaxis.limits = [-60, 60]
            if pm == 'E':
                plot.xaxis.limits = [E0-dE,E0+dE]
            if pm == 'Zp':
                plot.xaxis.limits = [-800, 800]

        if od == 'S3':
            if pm == 'EZp':
                plot.yaxis.limits = [E0-2*dE,E0+2*dE]
                plot.xaxis.limits = [-1200, 1200]
            else:
                plot.yaxis.limits = [-3.5, 3.5]
            if pm == 'X':
               plot.xaxis.limits = [-7, 7]
            if pm == 'E':
                plot.xaxis.limits = [E0-0.1,E0+0.1]
            if pm == 'Zp':
                plot.xaxis.limits = [-1200, 1200]

        if od == 'S4':
            if pm == 'EZp':
                plot.yaxis.limits = [E0-2*dE,E0+2*dE]
                #plot.xaxis.limits = [-20000, 20000]
            if pm == 'X':
                pass
#               plot.xaxis.limits = [-3.5, 3.5]
            if pm == 'E':
                plot.xaxis.limits = [E0-2*dE,E0+2*dE]
                plot.xaxis.offset = E0
#            if pm == 'Zp':
#                plot.xaxis.limits = [-20000, 20000]

        if od in ['S0','S1']:
            plot.caxis.limits = [E0-dE,E0+dE]
        if od in ['S2']:
            plot.caxis.limits = [E0-2*dE,E0+2*dE]
        if od in ['S3','S4']:
            ftr = 1
            #plot.caxis.limits = [(E0-dE)*ftr,(E0+dE)*ftr]
            plot.caxis.limits = [E0-dE,E0+dE]
            plot.caxis.factor = ftr
            plot.caxis.unit = r"eV"
            plot.caxis.fwhmFormatStr = '%1.4f'
            plot.caxis.label = r"energy"
        plot.caxis.offset = E0
        plot.xaxis.fwhmFormatStr = '%1.6f'
        plot.yaxis.fwhmFormatStr = '%1.6f'
#        plot.caxis.fwhmFormatStr = '%1.6f'

#    Plot02 = xrtplot.XYCPlot(
#        beam=r"screen01beamLocal01",
#        xaxis=xrtplot.XYCAxis(
#            label=r"x",
#            fwhmFormatStr=r"%.3f"),
#        yaxis=xrtplot.XYCAxis(
#            label=r"z",
#            fwhmFormatStr=r"%.3f"),
#        caxis=xrtplot.XYCAxis(
#            label=r"energy",
#            unit=r"eV",
#            fwhmFormatStr=r"%.3f"),
#        title=r"Screen between the two crystal")

    return plots

# important!
raycing.run.run_process = run_process



if __name__ == '__main__':

    try:
        step = float(sys.argv[1])
        det_ang_offset = radians(step)
    except:
        step = None

    try:
        step2 = float(sys.argv[2])
        Estep = step2
    except:
        step2 = None

    BL = build_beamline()
    BL.alignE = E0
    if glow3D:
        BL.glow(scale=[10, 1, 1],centerAt='Screen5')
        sys.exit()

    stats = beamline_stats(BL)
    # print(stats)

    if just_calc:
        sys.exit()

    plots = define_plots(ss, mode, step, step2)
    xrtrun.run_ray_tracing(beamLine=BL, repeats=repeats, plots=plots)

#    if step:
#        fout = f'{dat_dir}{title}_{step:.3f}'
#        if step2 is not None:
#            fout += f'_{step2:.3f}'
#        for plot in plots:
#            sn = plot.saveName[:-4].split('/')[-1]
#            od = sn.split('_')[2]
#            pm = sn.split('_')[1]
#            if od == 'S4' and pm == 'D':
#                y = plot.yaxis.total1D
#                x = plot.yaxis.binEdges[:-1]
#                np.savetxt(fout+'.dat', np.array([x,y]).T, header=stats)

    fout = f'{dat_dir}{title}_{abs(degrees(det_ang_offset)):.1f}'
    fout += f'_{enstep:.3f}'
    for plot in plots:
        sn = plot.saveName[:-4].split('/')[-1]
        od = sn.split('_')[2]
        pm = sn.split('_')[1]
        if od == 'S4' and pm == 'D':
            y = plot.yaxis.total1D
            x = plot.yaxis.binEdges[:-1]
            np.savetxt(fout+'.dat', np.array([x,y]).T, header=stats)
            print('\nData saved at {}.dat\n'.format(fout))