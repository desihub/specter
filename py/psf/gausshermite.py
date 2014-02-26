#!/usr/bin/env python
"""
GaussHermitePSF - PSF modeled with 2D Gauss-Hermite polynomials

Stephen Bailey
December 2013
"""

import sys
import os
import numpy as N
from scipy import special as sp
from numpy.polynomial.legendre import legval2d, legval

import fitsio
from specter.psf import PSF
from specter.util import TraceSet

def pgh(x, m=0, xc=0., sigma=1., dxlo=-0.5, dxhi=0.5):
    """
    Pixel-integrated (probabilist) Gauss-Hermite function.

    Arguments:
      x: pixel-center baseline array
      m: order of Hermite polynomial multiplying Gaussian core
      xc: sub-pixel position of Gaussian centroid relative to x baseline
      sigma: sigma parameter of Gaussian core in units of pixels
      dxlo, dxhi: pixel boundary low and high offsets, defauting to -/+0.5

    Written: Adam S. Bolton, U. of Utah, fall 2010
    """
    u = (x - xc) / sigma
    ulo = (x + dxlo - xc) / sigma
    uhi = (x + dxhi - xc) / sigma
    if m > 0:
        hiterm = - sp.hermitenorm(m-1)(uhi) * N.exp(-0.5 * uhi**2) / N.sqrt(2. * N.pi)
        loterm = - sp.hermitenorm(m-1)(ulo) * N.exp(-0.5 * ulo**2) / N.sqrt(2. * N.pi)
        # poly = sp.hermitenorm(m-1).coeffs
        # hiterm = - N.polyval(poly, uhi) * N.exp(-0.5 * uhi**2) / N.sqrt(2. * N.pi)
        # loterm = - N.polyval(poly, ulo) * N.exp(-0.5 * ulo**2) / N.sqrt(2. * N.pi)
        return hiterm - loterm
    else:
        return 0.5 * (sp.erf(uhi/N.sqrt(2.)) - sp.erf(ulo/N.sqrt(2.)))

class GaussHermitePSF(PSF):
    """
    Model PSF with a linear interpolation of high resolution sampled spots
    """
    def __init__(self, filename):
        """
        Initialize SpotGridPSF from input file
        
        See specter.psf.PSF for futher details
        """
        #- Use PSF class to Load Generic PSF values (x, y, wavelength, ...)
        ### PSF.__init__(self, filename)
        
        #- Check that this file is a recent Gauss Hermite PSF
        fx = fitsio.FITS(filename)
        hdr = fx[1].read_header()
        if 'PSFTYPE' not in hdr:
            raise ValueError, 'Missing PSFTYPE keyword'
            
        if hdr['PSFTYPE'] != 'GAUSS-HERMITE2':
            raise ValueError, 'PSFTYPE %s is not GAUSS-HERMITE2' % hdr['PSFTYPE']
            
        if 'PSFVER' not in hdr:
            raise ValueError, "PSFVER missing; this version not supported"
            
        # if hdr['PSFVER'] < '1':
        #     raise ValueError, "Only GAUSS-HERMITE versions 1.0 and greater are supported"
            
        #- HACK
        self.nspec = hdr['FIBERMAX'] - hdr['FIBERMIN'] + 1
        
        #- Other necessary keywords
        self.npix_x = hdr['NPIX_X']
        self.npix_y = hdr['NPIX_Y']
        
        #- Load the parameters
        #- Hack: need to assemble bundles
        hdr = fx[1].read_header()
        data = fx[1].read()

        self._coeff = dict()
        self._polyparams = hdr
        for p in data:
            domain = (p['WAVEMIN'], p['WAVEMAX'])
            for p in data:
                name = p['PARAM'].strip()
                self._coeff[name] = TraceSet(p['COEFF'], domain=domain)
        
        #- Pull out x and y as special tracesets
        self._x = self._coeff['X']
        self._y = self._coeff['Y']

        #- Create inverse y -> wavelength mapping
        self._w = self._y.invert()
        self._wmin = N.min(self.wavelength(None, 0))
        self._wmax = N.max(self.wavelength(None, self.npix_y-1))
                
        #- Filled only if needed
        self._xsigma = None
        self._ysigma = None
        
    def _xypix(self, ispec, wavelength):
                
        # x, y = self.xy(ispec, wavelength)
        x = self._coeff['X'].eval(ispec, wavelength)
        y = self._coeff['Y'].eval(ispec, wavelength)
        
        #- CCD pixel ranges
        hsizex = self._polyparams['HSIZEX']
        hsizey = self._polyparams['HSIZEY']
        xccd = N.arange(int(x-hsizex), int(x+hsizex))
        yccd = N.arange(int(y-hsizey), int(y+hsizey))
        dx = xccd - x
        dy = yccd - y
        
        #- Extract coefficients
        degx1 = self._polyparams['GHDEGX']
        degy1 = self._polyparams['GHDEGY']
        degx2 = self._polyparams['GHDEGX2']
        degy2 = self._polyparams['GHDEGY2']
        sigx1 = self._coeff['GHSIGX'].eval(ispec, wavelength)
        sigx2 = self._coeff['GHSIGX2'].eval(ispec, wavelength)
        sigy1 = self._coeff['GHSIGY'].eval(ispec, wavelength)
        sigy2 = self._coeff['GHSIGY2'].eval(ispec, wavelength)        
        
        #- Background tail image
        tailxsca = self._coeff['TAILXSCA'].eval(ispec, wavelength)
        tailysca = self._coeff['TAILYSCA'].eval(ispec, wavelength)
        tailamp = self._coeff['TAILAMP'].eval(ispec, wavelength)
        tailcore = self._coeff['TAILCORE'].eval(ispec, wavelength)
        tailinde = self._coeff['TAILINDE'].eval(ispec, wavelength)

        #- Make tail image
        #- TODO: This could be done faster by removing loops
        img = N.zeros((len(yccd), len(xccd)))
        for i, dyy in enumerate(dy):
            for j, dxx in enumerate(dx):
                r2 = (dxx*tailxsca)**2 + (dyy*tailysca)**2
                img[i,j] = tailamp * r2 / (tailcore**2 + r2)**(1+tailinde/2.0)

        #- Create 1D GaussHermite functions in x and y
        xfunc1 = [pgh(xccd, i, x, sigma=sigx1) for i in range(degx1+1)]
        yfunc1 = [pgh(yccd, i, y, sigma=sigy1) for i in range(degy1+1)]        
        for i in range(degx1+1):
            for j in range(degy1+1):
                spot1 = N.outer(yfunc1[j], xfunc1[i])
                c1 = self._coeff['GH-{}-{}'.format(i,j)].eval(ispec, wavelength)
                img += c1 * spot1

        xfunc2 = [pgh(xccd, i, x, sigma=sigx2) for i in range(degx2+1)]
        yfunc2 = [pgh(yccd, i, y, sigma=sigy2) for i in range(degy2+1)]
        img2 = N.zeros((len(yccd), len(xccd)))
        for i in range(degx2+1):
            for j in range(degy2+1):
                spot2 = N.outer(yfunc2[j], xfunc2[i])
                c2 = self._coeff['GH2-{}-{}'.format(i,j)].eval(ispec, wavelength)
                img += c2 * spot2

        xslice = slice(xccd[0], xccd[-1]+1)
        yslice = slice(yccd[0], yccd[-1]+1)
        return xslice, yslice, img
        
#-------------------------------------------------------------------------
# for igx in range(p['GHDEGX']+1):
#     for igy in range(p['GHDEGY']+1):
#         subplot(5,5,igx*5 + igy+1)
#         ghcoeff = legval2d(xx, yy, c[:,:,igx,igy])
#         print igx, igy, ghcoeff
#         # imshow( ghcoeff * N.outer(yfuncs[igy], xfuncs[igx]), vmin=0, vmax=0.5 )
#         img = N.outer(yfuncs[igy], xfuncs[igx])
#         vrange = N.max(N.abs(img))
#         imshow( img, vmin=-vrange, vmax=vrange)

