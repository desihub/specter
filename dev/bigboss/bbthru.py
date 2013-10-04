#!/usr/bin/env python

"""
Collect BigBOSS throughput pieces from bbspecsim into
Specter throughput format.

TODO: Factor our parameters at the top into a configuration file.

Stephen Bailey, LBL
January 2013
"""

import sys
import os
import numpy as N
import fitsio

from scipy.interpolate import InterpolatedUnivariateSpline as SplineFit

#--- Input Parameters ---
#- TODO: Move these to a configuration file

#- Wavelength range and resolution
wmin, wmax, dw = 3500, 9900, 0.1

#- Fiber parameters
fiber_size_um  = 120.0          #- diameter in microns
fiber_frd_loss = 0.970          #- Focal Ratio Degradation
fiber_connection_loss = 0.975   #- Loss from fiber connections
fiberlen = 40.0                 #- meters

#- Default exposure time
exptime = 15*60.0  #- 15 minutes -> seconds

#- Input directories, corrector and spectrograph versions
bbdir = os.environ['BBSPECSIM_DIR']
bbcorr = 'BB_CORR_20120402a'
bbspec = 'BB_SPEC_20120428difdet'

extinction_file = bbdir + '/sky/ZenExtinct-KPNO-FTS.fits'
corrector_file = '%s/designs/%s/%s.fits' % (bbdir, bbcorr, bbcorr)
#- use "design_file % n" to get file for spectrograph n=1,2,3
design_file = '%s/designs/%s/%s_%%d.fits.gz' % (bbdir, bbspec, bbspec)

#--- Data to fill ---
thru = dict()
ww = N.arange(wmin, wmax+dw/2, dw)
nw = len(ww)
thru['wavelength'] = ww

#- Vectors of throughputs common to all spectrographs
thru['extinction'] = None    #- Atmosphere
thru['optics']     = None    #- Telescope + corrector but not spectrographs
thru['fiberinput'] = None    #- Fiber input geometric losses
thru['fiber']      = None    #- Fiber optics

#- Individual spectrographs
thru['B'] = N.zeros(nw)
thru['R'] = N.zeros(nw)
thru['I'] = N.zeros(nw)

#- Scalars
thru['effarea'] = None
thru['fiberdia'] = None
thru['exptime'] = exptime

#--- Load throughputs from bbspecsim files ---

#- Atmospheric extinction: magnitudes per airmass
#- throughput = 10**(-0.4*airmass*extinction)
atm = fitsio.read(extinction_file, 1, lower=True)
s = SplineFit(atm['wavelength'], atm['extinction'])
thru['extinction'] = s(ww)

#- Telescope + Corrector model
tele = fitsio.read(corrector_file, upper=True)
tw = tele['WAVE'][0]
tx = tele['MIRROREF'][0] * tele['CORRTRANS'][0] * tele['ARTRANS'][0]**tele['NSURF'][0]
s = SplineFit(tw, tx)
thru['optics'] = s(ww)
thru['effarea'] = N.pi * (tele['MIRRORDIAM'][0]*100)**2 / 4 * tele['VIGN'][0]  #- cm^2
thru['fiberdia'] = fiber_size_um / tele['FPSCALE'][0]  #- microns -> arcsec

#- Fiber input model
#- WARNING: this is seeing and object shape dependent, and uses a
#- file leftover from a previous bbspecsim run.
fin = N.loadtxt(bbdir+'/throughput/fiberee.txt').T
s = SplineFit(fin[0], fin[1])
thru['fiberinput'] = s(ww)

#- Fiber attenuation (dB/km) = 10*log10(in/out)
fatt = N.loadtxt(bbdir+'/throughput/FiberAttenuation.txt').T
fw = fatt[0]*10
fa = fatt[1]
fx = 10**(-fa/10.0 * (fiberlen/1000.0) )
thru['fiber'] = N.interp(ww, fw, fx)  #- linear interp instead of ringy spline
thru['fiber'] *= fiber_frd_loss
thru['fiber'] *= fiber_connection_loss

#- Spectrograph efficiencies
specname = {1:'B', 2:'R', 3:'I'}
for ispec in (1,2,3):
    spec = fitsio.read(design_file % ispec, 2, lower=True)
    spw = spec['lambda'][0] * 10  #- nm -> A
    s = SplineFit(spw, spec['spectotal'][0]*spec['detector'][0])
    #- Don't extrapolate spline
    ii = (spw[0] <= ww) & (ww <= spw[-1])
    thru[specname[ispec]][ii] = s(ww[ii])

#--- Write one fits file per spectrograph ---
hdr = list()
hdr.append(dict(name='EXTNAME', value='THROUGHPUT'))
hdr.append(dict(name='EXPTIME', value=thru['exptime'],
    comment='Default exposure time [seconds]'))
hdr.append(dict(name='EFFAREA', value=thru['effarea'],
    comment='effective mirror area [cm^2]'))
hdr.append(dict(name='FIBERDIA', value=thru['fiberdia'],
    comment='Fiber diameter [arcsec]'))

wavelength = thru['wavelength']
extinction = thru['extinction']
fiberinput = thru['fiberinput']
for camera in specname.values():
    throughput = thru['optics'] * thru['fiber'] * thru[camera]
    thru[camera+'tot'] = throughput * thru['fiberinput']
    ii = N.where(throughput > 0)
    data = N.rec.fromarrays((wavelength[ii], extinction[ii], throughput[ii], fiberinput[ii]),
                            names='wavelength,extinction,throughput,fiberinput')
    outfile = 'bbthru-%s.fits' % camera
    fitsio.write(outfile, data, header=hdr, clobber=True)

#--- DEBUG ---
import pylab as P
P.ion()
P.rcParams['legend.fontsize'] = 12
P.plot(ww, thru['optics'], 'k-', lw=1, alpha=0.5, label='Tele+Corr Optics')
P.plot(ww, thru['fiberinput'], 'g--', lw=2, label='Fiber Input (%.2f arcsec)' % thru['fiberdia'])
P.plot(ww, thru['fiber'], 'm:', lw=3, label='Fibers (%.1f m)' % fiberlen)
P.fill_between(ww, thru['Btot'], color='b', alpha=0.5)
P.plot(ww, thru['B'], color='b', label='B spectro total')
P.fill_between(ww, thru['Rtot'], color='r', alpha=0.5)
P.plot(ww, thru['R'], color='r', label='R spectro total')
P.fill_between(ww, thru['Itot'], color='k', alpha=0.5)
P.plot(ww, thru['I'], color='k', label='I spectro total')
P.legend(loc='upper left', frameon=False)
P.title('BigBOSS Throughputs (no atmosphere)')
P.xlabel('Wavelength [A]')
P.ylabel('Throughput')
P.savefig('bbthru.png', dpi=80)
import IPython
IPython.embed()
#--- DEBUG ---
