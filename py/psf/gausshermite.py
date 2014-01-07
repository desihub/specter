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
        PSF.__init__(self, filename)
        
        #- Check that this file is a recent Gauss Hermite PSF
        fx = fitsio.FITS(filename)
        hdr = fx[0].read_header()
        if 'PSFTYPE' not in hdr:
            raise ValueError, 'Missing PSFTYPE keyword'
            
        if hdr['PSFTYPE'] != 'GAUSS-HERMITE':
            raise ValueError, 'PSFTYPE %s is not GAUSS-HERMITE' % hdr['PSFTYPE']
            
        if 'PSFVER' not in hdr:
            raise ValueError, "PSFVER missing; this version not supported"
            
        if hdr['PSFVER'] < '2':
            raise ValueError, "Only GAUSS-HERMITE versions 2.0 and greater are supported"
            
        nspec = hdr['NSPEC']
        
        #- Load the parameters
        self._coeff = dict()
        self._polyparams = dict()
        for bundle in fx[3:]:
            hdr = bundle.read_header()
            data = bundle.read()
            ### for ispec in range(hdr['FIBERMIN'], hdr['FIBERMAX']+1):
            for ispec in range(nspec):           #- HACK!
                self._coeff[ispec] = data
                self._polyparams[ispec] = hdr
        
    def _xypix(self, ispec, wavelength):
        c = self._coeff[ispec]
        p = self._polyparams[ispec]
        x, y = self.xy(ispec, wavelength)
        
        #- Convert x,y to range [-1,1] for Legendre
        xx = 2.0*(x - p['LXMIN']) / (p['LXMAX'] - p['LXMIN']) - 1.0
        yy = 2.0*(y - p['LYMIN']) / (p['LYMAX'] - p['LYMIN']) - 1.0
        
        #- CCD pixel ranges
        halfsize = 5
        xccd = N.arange(int(x-halfsize), int(x+halfsize))
        yccd = N.arange(int(y-halfsize), int(y+halfsize))
        dx = xccd - x
        dy = yccd - y
        
        #- Build Gauss-Hermite functions
        xfuncs = N.zeros((p['GHDEGX']+1, len(xccd)), dtype=float)
        yfuncs = N.zeros((p['GHDEGY']+1, len(yccd)), dtype=float)
        for i in range(p['GHDEGX']+1):
            xfuncs[i] = pgh(xccd, i, x, p['GHSIGX'])
        for i in range(p['GHDEGY']+1):
            yfuncs[i] = pgh(yccd, i, y, p['GHSIGY'])
        
        # psfimg = N.zeros((len(yccd), len(xccd)))
        
        psfimg = list()
        ghcoeff = list()
        for igx in range(p['GHDEGX']+1):
            for igy in range(p['GHDEGY']+1):
                ghcoeff.append( legval2d(yy, xx, c[igy,igx]) )
                psfimg.append( N.outer(yfuncs[igy], xfuncs[igx]) )
        
        psfimg = N.array(psfimg)
        ghcoeff = N.array(ghcoeff)
        
        img = psfimg.T.dot(ghcoeff).T
        
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

