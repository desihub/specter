#!/usr/bin/env python
"""
GaussHermitePSF - PSF modeled with 2D Gauss-Hermite polynomials as generated
by the specex package at https://github.com/julienguy/specex

Stephen Bailey
December 2013
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
from scipy import special as sp

from astropy.io import fits
from specter.psf import PSF
from specter.util import TraceSet, outer, legval_numba, custom_hermitenorm, custom_erf

import numba

class GaussHermitePSF(PSF):
    """
    Model PSF with two central Gauss-Hermite cores with different sigmas
    plus power law wings.
    """
    def __init__(self, filename):
        """
        Initialize GaussHermitePSF from input file
        """        
        #- Check that this file is a current generation Gauss Hermite PSF
        fx = fits.open(filename, memmap=False)

        #- Read primary header
        phdr = fx[0].header
        if 'PSFTYPE' not in phdr:
            raise ValueError('Missing PSFTYPE keyword')
        if phdr['PSFTYPE'] != 'GAUSS-HERMITE':
            raise ValueError('PSFTYPE {} is not GAUSS-HERMITE'.format(phdr['PSFTYPE']))
        if 'PSFVER' not in phdr:
            raise ValueError("PSFVER missing; this version not supported")
        PSFVER = float(phdr["PSFVER"])
        
        if PSFVER<3 :
            psf_hdu = 1
        else :
            psf_hdu = "PSF"
        
        self._polyparams = hdr = dict(fx[psf_hdu].header)
        if 'PSFTYPE' not in hdr:
            raise ValueError('Missing PSFTYPE keyword')
            
        if hdr['PSFTYPE'] != 'GAUSS-HERMITE':
            raise ValueError('PSFTYPE {} is not GAUSS-HERMITE'.format(hdr['PSFTYPE']))
            
        if 'PSFVER' not in hdr:
            raise ValueError("PSFVER missing; this version not supported")
            
        if hdr['PSFVER'] < '1':
            raise ValueError("Only GAUSS-HERMITE versions 1.0 and greater are supported")
            
        #- Calculate number of spectra from FIBERMIN and FIBERMAX (inclusive)
        self.nspec = hdr['FIBERMAX'] - hdr['FIBERMIN'] + 1
        
        #- Other necessary keywords
        self.npix_x = hdr['NPIX_X']
        self.npix_y = hdr['NPIX_Y']

        #- PSF model error
        if 'PSFERR' in hdr:
            self.psferr = hdr['PSFERR']
        else:
            self.psferr = 0.01

        #- Load the parameters into self.coeff dictionary keyed by PARAM
        #- with values as TraceSets for evaluating the Legendre coefficients
        data = fx[psf_hdu].data
        self.coeff = dict()
        if PSFVER<3 : # old format
            
            for p in data:
                domain = (p['WAVEMIN'], p['WAVEMAX'])
                for p in data:
                    name = p['PARAM'].strip()
                    self.coeff[name] = TraceSet(p['COEFF'], domain=domain)
    
            #- Pull out x and y as special tracesets
            self._x = self.coeff['X']
            self._y = self.coeff['Y']

        else : # new format
            
            domain = (hdr['WAVEMIN'], hdr['WAVEMAX'])
            for p in data:
                name = p['PARAM'].strip()
                self.coeff[name] = TraceSet(p['COEFF'], domain=domain)
            
            self._x = TraceSet(fx["XTRACE"].data, domain=(fx['XTRACE'].header["WAVEMIN"], fx['XTRACE'].header['WAVEMAX']))
            self._y = TraceSet(fx["YTRACE"].data, domain=(fx['YTRACE'].header["WAVEMIN"], fx['YTRACE'].header['WAVEMAX']))
            
        
        #- Create inverse y -> wavelength mapping
        self._w = self._y.invert()

        #- Cache min/max wavelength per fiber at pixel edges
        self._wmin_spec = self.wavelength(None, -0.5)
        self._wmax_spec = self.wavelength(None, self.npix_y-0.5)
        self._wmin = np.min(self._wmin_spec)
        self._wmin_all = np.max(self._wmin_spec)
        self._wmax = np.max(self._wmax_spec)
        self._wmax_all = np.min(self._wmax_spec)

        #- Filled only if needed
        self._xsigma = None
        self._ysigma = None

        #create dict to hold legval cached data
        self.legval_dict = None

        #- Cache hermitenorm polynomials so we don't have to create them
        #- every time xypix is called
        self._hermitenorm = list()
        maxdeg = max(hdr['GHDEGX'], hdr['GHDEGY'])
        for i in range(maxdeg+1):
            self._hermitenorm.append( sp.hermitenorm(i) )

        fx.close()

    def _xypix(self, ispec, wavelength, ispec_cache=None, iwave_cache=None):
        """
        Two branches of this function which does legendre series fitting
        First branch = no caching, will recompute legval every time eval is used (slow)
        Second branch = yes caching, will look up precomputed values of legval (faster)

        Arguments:
            ispec: the index of each spectrum
            wavelength: the wavelength at which to evaluate the legendre series
            
        Optional arguments:
            ispec_cache: the index of each spectrum which starts again at 0 for each patch
            iwave_cache: the index of each wavelength which starts again at 0 for each patch
        """

        #trying to avoid duplicating code between branches, will split into
        #cached branch and non-cached branch depending on self.legval_dict

        #we need this for the original indexing so skip looking up cached values here
        x = self._x.eval(ispec, wavelength)
        y = self._y.eval(ispec, wavelength)

        #- CCD pixel ranges
        hsizex = self._polyparams['HSIZEX']
        hsizey = self._polyparams['HSIZEY']
        xccd = np.arange(int(x-hsizex+0.5), int(x+hsizex+1.5))
        yccd = np.arange(int(y-hsizey+0.5), int(y+hsizey+1.5))

        dx = xccd - x
        dy = yccd - y
        nx = len(dx)
        ny = len(dy)

        degx1 = self._polyparams['GHDEGX']
        degy1 = self._polyparams['GHDEGY']

        if self.legval_dict is None:
            #non-cached branch
            #print("self.legval_dict is None, hitting non-cached branch")

            #- Extract GH degree and sigma coefficients for convenience
            sigx1 = self.coeff['GHSIGX'].eval(ispec, wavelength)
            sigy1 = self.coeff['GHSIGY'].eval(ispec, wavelength)
        
            #- Background tail image
            tailxsca = self.coeff['TAILXSCA'].eval(ispec, wavelength)
            tailysca = self.coeff['TAILYSCA'].eval(ispec, wavelength)
            tailamp = self.coeff['TAILAMP'].eval(ispec, wavelength)
            tailcore = self.coeff['TAILCORE'].eval(ispec, wavelength)
            tailinde = self.coeff['TAILINDE'].eval(ispec, wavelength)
 
        else:
            #cached branch
            #print("self.legval_dict is not None, hitting cached branch")

            #- Extract GH degree and sigma coefficients for convenience
            sigx1 = self.legval_dict['GHSIGX'][ispec_cache, iwave_cache]
            sigy1 = self.legval_dict['GHSIGY'][ispec_cache, iwave_cache]

            #- Background tail image
            tailxsca = self.legval_dict['TAILXSCA'][ispec_cache, iwave_cache]
            tailysca = self.legval_dict['TAILYSCA'][ispec_cache, iwave_cache]
            tailamp = self.legval_dict['TAILAMP'][ispec_cache, iwave_cache]
            tailcore = self.legval_dict['TAILCORE'][ispec_cache, iwave_cache]
            tailinde = self.legval_dict['TAILINDE'][ispec_cache, iwave_cache]

        
        #- Make tail image (slow version)
        # img = np.zeros((len(yccd), len(xccd)))
        # for i, dyy in enumerate(dy):
        #     for j, dxx in enumerate(dx):
        #         r2 = (dxx*tailxsca)**2 + (dyy*tailysca)**2
        #         img[i,j] = tailamp * r2 / (tailcore**2 + r2)**(1+tailinde/2.0)

        #- Make tail image (faster, less readable version)
        #- r2 = normalized distance from center of each pixel to PSF center
        r2 = np.tile((dx*tailxsca)**2, ny).reshape(ny, nx) + \
             np.repeat((dy*tailysca)**2, nx).reshape(ny, nx)
        tails = tailamp*r2 / (tailcore**2 + r2)**(1+tailinde/2.0)
        
        #- Create 1D GaussHermite functions in x and y
        xfunc1 = [pgh(xccd, i, x, sigma=sigx1) for i in range(degx1+1)]
        yfunc1 = [pgh(yccd, i, y, sigma=sigy1) for i in range(degy1+1)]        
       
        #- Create core PSF image
        core1 = np.zeros((ny, nx))
        spot1 = np.empty_like(core1)

        if self.legval_dict is None:
            #non-cached branch

            for i in range(degx1+1):
                for j in range(degy1+1):
                    c1 = self.coeff['GH-{}-{}'.format(i,j)].eval(ispec, wavelength)
                    outer(yfunc1[j], xfunc1[i], out=spot1)
                    core1 += c1 * spot1

        else:
            #cached branch

            #for i in range(degx1+1):
            #    for j in range(degy1+1):
            #        #see if we can squeeze out some extra speed by looking up the values
            #        c1 = self.legval_dict['GH-{}-{}'.format(i,j)][ispec_cache, iwave_cache]
            #        outer(yfunc1[j], xfunc1[i], out=spot1)
            #        core1 += c1 * spot1
            #try our own function
            #but first some prep
            npfx = np.array(xfunc1)
            npfy = np.array(yfunc1)

            #store values of c1
            c1_array = np.zeros((ny, nx))

            #get legval_value
            for i in range(degx1+1):
                for j in range(degy1+1):
                    #see if we can squeeze out some extra speed by looking up the values
                    c1_array[i,j] = self.legval_dict[self.core_keys[i][j]][ispec_cache, iwave_cache]

            core1 = generate_core(degx1, degy1, npfx, npfy, spot1, core1, c1_array)

        #- Zero out elements in the core beyond 3 sigma
        #- Only for GaussHermite2
        # ghnsig = self.coeff['GHNSIG'].eval(ispec, wavelength)
        # r2 = np.tile((dx/sigx1)**2, ny).reshape(ny, nx) + \
        #      np.repeat((dy/sigy1)**2, nx).reshape(ny, nx)
        
        # core1 *= (r2<ghnsig**2)
        
        #this code will not work, needs to be modified for new pgh!
        #- Add second wider core Gauss-Hermite term        
        # xfunc2 = [self._pgh(xccd, i, x, sigma=sigx2) for i in range(degx2+1)]
        # yfunc2 = [self._pgh(yccd, i, y, sigma=sigy2) for i in range(degy2+1)]
        # core2 = np.zeros((ny, nx))
        # for i in range(degx2+1):
        #     for j in range(degy2+1):
        #         spot2 = outer(yfunc2[j], xfunc2[i])
        #         c2 = self.coeff['GH2-{}-{}'.format(i,j)].eval(ispec, wavelength)
        #         core2 += c2 * spot2

        #- Clip negative values and normalize to 1.0
        img = core1 + tails
        ### img = core1 + core2 + tails

        img = img.clip(0.0)
        img /= np.sum(img)

        xslice = slice(xccd[0], xccd[-1]+1)
        yslice = slice(yccd[0], yccd[-1]+1)
        return xslice, yslice, img
        # return xslice, yslice, (core1, core2, tails)
        
    def _gh(self, x, m=0, xc=0.0, sigma=1.0):
        """
        return Gauss-Hermite function value, NOT integrated, for display of PSF.
        Arguments:
          x: coordinates baseline array
          m: order of Hermite polynomial multiplying Gaussian core
          xc: sub-pixel position of Gaussian centroid relative to x baseline
          sigma: sigma parameter of Gaussian core in units of pixels
        
        """

        
        u = (x-xc) / sigma
        
        if m > 0:
            return self._hermitenorm[m](u) * np.exp(-0.5 * u**2) / np.sqrt(2. * np.pi)
        else:                        
            return np.exp(-0.5 * u**2) / np.sqrt(2. * np.pi)       


    def xsigma(self, ispec, wavelength):
        """
        Return Gaussian sigma of PSF spot in cross-dispersion direction
        in CCD pixel units.
        
        ispec : spectrum index
        wavelength : scalar or vector wavelength(s) to evaluate spot sigmas
        """
        return self.coeff['GHSIGX'].eval(ispec, wavelength)
    
    def ysigma(self, ispec, wavelength):
        """
        Return Gaussian sigma of PSF spot in dispersion direction
        in CCD pixel units.
        
        ispec : spectrum index
        wavelength : scalar or vector wavelength(s) to evaluate spot sigmas
        """
        return self.coeff['GHSIGY'].eval(ispec, wavelength)
    

    def _value(self,x,y,ispec, wavelength):
        
        """
        return PSF value (same shape as x and y), NOT integrated, for display of PSF.

        Arguments:
          x: x-coordinates baseline array
          y: y-coordinates baseline array (same shape as x)
          ispec: fiber
          wavelength: wavelength
       
        """

        # x, y = self.xy(ispec, wavelength)
        xc = self._x.eval(ispec, wavelength)
        yc = self._y.eval(ispec, wavelength)
                
        #- Extract GH degree and sigma coefficients for convenience
        degx1 = self._polyparams['GHDEGX']
        degy1 = self._polyparams['GHDEGY']
        sigx1 = self.coeff['GHSIGX'].eval(ispec, wavelength)
        sigy1 = self.coeff['GHSIGY'].eval(ispec, wavelength)
        
        #- Background tail image
        tailxsca = self.coeff['TAILXSCA'].eval(ispec, wavelength)
        tailysca = self.coeff['TAILYSCA'].eval(ispec, wavelength)
        tailamp = self.coeff['TAILAMP'].eval(ispec, wavelength)
        tailcore = self.coeff['TAILCORE'].eval(ispec, wavelength)
        tailinde = self.coeff['TAILINDE'].eval(ispec, wavelength)
        
        
        #- Make tail image (faster, less readable version)
        r2 = ((x-xc)*tailxsca)**2+((y-yc)*tailysca)**2
        tails = tailamp*r2 / (tailcore**2 + r2)**(1+tailinde/2.0)
        
        #- Create 1D GaussHermite functions in x and y
        xfunc1 = [self._gh(x, i, xc, sigma=sigx1) for i in range(degx1+1)]
        yfunc1 = [self._gh(y, i, yc, sigma=sigy1) for i in range(degy1+1)]        
        
        
        #- Create core PSF image
        core1 = np.zeros(x.shape)
        for i in range(degx1+1):
            for j in range(degy1+1):
                c1 = self.coeff['GH-{}-{}'.format(i,j)].eval(ispec, wavelength)
                core1 += c1 * yfunc1[j]*xfunc1[i]
        
        #- Clip negative values and normalize to 1.0
        img = core1 + tails
        ### img = core1 + core2 + tails

        img = img.clip(0.0)
        img /= np.sum(img)

        return img

    def cache_params(self, specrange, wavelengths):
        #store in a dict
        self.legval_dict = dict()
        #store keys in a list
        self.core_keys = list()
        for key in self.coeff.keys():
            self.legval_dict[key] = self.coeff[key].eval(specrange, wavelengths)
        #some extra steps to cache what we need for the core PSF image
        degx1 = self._polyparams['GHDEGX']
        degy1 = self._polyparams['GHDEGY']
        for i in range(degx1+1):
            self.core_keys.append(list())
            for j in range(degy1+1):
                self.core_keys[-1].append('GH-{}-{}'.format(i,j))
                self.legval_dict['GH-{}-{}'.format(i,j)]=self.coeff[self.core_keys[i][j]].eval(specrange, wavelengths)

@numba.jit(nopython=True, cache=False)
def pgh(x, m=0, xc=0.0, sigma=1.0):
    """
    Pixel-integrated (probabilist) Gauss-Hermite function.
    Arguments:
      x: pixel-center baseline array
      m: order of Hermite polynomial multiplying Gaussian core
      xc: sub-pixel position of Gaussian centroid relative to x baseline
      sigma: sigma parameter of Gaussian core in units of pixels
    Uses the relationship
    Integral{ H_k(x) exp(-0.5 x^2) dx} = -H_{k-1}(x) exp(-0.5 x^2) + const
    Written: Adam S. Bolton, U. of Utah, fall 2010
    Adapted for efficiency by S. Bailey while dropping generality
    modified from the orig _pgh to enable jit-compiling
    --> no longer passing in self, calling custom numba functions in util
    """

    #- Evaluate H[m-1] at half-pixel offsets above and below x
    dx = x-xc-0.5
    u = np.concatenate( (dx, dx[-1:]+0.5) ) / sigma
        
    if m > 0:
        y  = -custom_hermitenorm(m-1,u) * np.exp(-0.5 * u**2) / np.sqrt(2. * np.pi)
        return (y[1:] - y[0:-1])
    else:            
        y = custom_erf(u/np.sqrt(2.))
        return 0.5 * (y[1:] - y[0:-1])

#new function to numba-ize the expensive part of _xypix
@numba.jit(nopython=True, cache=False)
def generate_core(degx1, degy1, npfx, npfy, spot1, core1, c1_array):
    for i in range(degx1+1):
        for j in range(degy1+1):
            outer(npfy[j], npfx[i], out=spot1)
            core1 += c1_array[i,j] * spot1
    return core1    
  


