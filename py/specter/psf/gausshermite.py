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
from specter.util import TraceSet, outer

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
        self._polyparams = hdr = dict(fx[1].header)
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
        data = fx[1].data
        self.coeff = dict()
        for p in data:
            domain = (p['WAVEMIN'], p['WAVEMAX'])
            for p in data:
                name = p['PARAM'].strip()
                self.coeff[name] = TraceSet(p['COEFF'], domain=domain)

        #- Pull out x and y as special tracesets
        self._x = self.coeff['X']
        self._y = self.coeff['Y']

        #- Create inverse y -> wavelength mapping
        self._w = self._y.invert()
        self._wmin = np.min(self.wavelength(None, 0))
        self._wmin_all = np.max(self.wavelength(None, 0))
        self._wmax = np.max(self.wavelength(None, self.npix_y-1))
        self._wmax_all = np.min(self.wavelength(None, self.npix_y-1))

        #- Filled only if needed
        self._xsigma = None
        self._ysigma = None

        #- Cache hermitenorm polynomials so we don't have to create them
        #- every time xypix is called
        self._hermitenorm = list()
        maxdeg = max(hdr['GHDEGX'], hdr['GHDEGY'])
        for i in range(maxdeg+1):
            self._hermitenorm.append( sp.hermitenorm(i) )

        #- No coefficients cached yet
        self._cache_ispec = None
        self._cache_hit = 0
        self._cache_miss = 0
        self._cache_coeff_keys = ['X', 'Y', 'GHSIGX', 'GHSIGY',
            'TAILXSCA', 'TAILYSCA', 'TAILAMP', 'TAILCORE', 'TAILINDE']
        self._cache_coeff = dict()
        degx1 = self._polyparams['GHDEGX']
        degy1 = self._polyparams['GHDEGY']
        for i in range(degx1+1):
            for j in range(degy1+1):
                key = 'GH-{}-{}'.format(i,j)
                self._cache_coeff_keys.append(key)

        fx.close()

    def _pgh(self, x, m=0, xc=0.0, sigma=1.0):
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
        """

        #- Evaluate H[m-1] at half-pixel offsets above and below x
        dx = x-xc-0.5
        u = np.concatenate( (dx, dx[-1:]+0.5) ) / sigma

        if m > 0:
            y = -self._hermitenorm[m-1](u) * np.exp(-0.5 * u**2) / np.sqrt(2. * np.pi)
            return (y[1:] - y[0:-1])
        else:
            y = sp.erf(u/np.sqrt(2.))
            return 0.5 * (y[1:] - y[0:-1])

    def cache_coeffs(self, ispecs, wavelengths):
        self._cache_wavelengths = wavelengths
        self._cache_coeff = dict()
        wavenorm = self.coeff['GH-0-0'].xnorm(wavelengths)
        for i in ispecs:
            self._cache_coeff[i] = dict()
            for key in self._cache_coeff_keys:
                self._cache_coeff[i][key] = self.coeff[key].eval(i, xnorm=wavenorm)

    def get_coeff(self, key, ispec, wave, iwave=None):
        if ispec in self._cache_coeff:
            if iwave is None:
                iwave = np.searchsorted(self._cache_wavelengths, wave)
            if iwave >= len(self._cache_wavelengths) or self._cache_wavelengths[iwave] != wave:
                self._cache_miss += 1
                return self.coeff[key].eval(ispec, wave)

            self._cache_hit += 1
            return self._cache_coeff[ispec][key][iwave]
        else:
            self._cache_miss += 1
            return self.coeff[key].eval(ispec, wave)

    def _xypix(self, ispec, wavelength):
        # x, y = self.xy(ispec, wavelength)
        wavenorm = self.coeff['GH-0-0'].xnorm(wavelength)
        # x = self.coeff['X'].eval(ispec, xnorm=wavenorm)
        # y = self.coeff['Y'].eval(ispec, xnorm=wavenorm)

        if hasattr(self, '_cache_wavelengths'):
            iwave = np.searchsorted(self._cache_wavelengths, wavelength)
        else:
            iwave = None

        x = self.get_coeff('X', ispec, wavelength, iwave)
        y = self.get_coeff('Y', ispec, wavelength, iwave)

        #- CCD pixel ranges
        hsizex = self._polyparams['HSIZEX']
        hsizey = self._polyparams['HSIZEY']
        xccd = np.arange(int(x-hsizex+0.5), int(x+hsizex+1.5))
        yccd = np.arange(int(y-hsizey+0.5), int(y+hsizey+1.5))
        dx = xccd - x
        dy = yccd - y
        nx = len(dx)
        ny = len(dy)

        #- Extract GH degree and sigma coefficients for convenience
        degx1 = self._polyparams['GHDEGX']
        degy1 = self._polyparams['GHDEGY']
        # sigx1 = self.coeff['GHSIGX'].eval(ispec, xnorm=wavenorm)
        # sigy1 = self.coeff['GHSIGY'].eval(ispec, xnorm=wavenorm)
        sigx1 = self.get_coeff('GHSIGX', ispec, wavelength, iwave)
        sigy1 = self.get_coeff('GHSIGY', ispec, wavelength, iwave)

        #- Background tail image
        # tailxsca = self.coeff['TAILXSCA'].eval(ispec, xnorm=wavenorm)
        # tailysca = self.coeff['TAILYSCA'].eval(ispec, xnorm=wavenorm)
        # tailamp = self.coeff['TAILAMP'].eval(ispec, xnorm=wavenorm)
        # tailcore = self.coeff['TAILCORE'].eval(ispec, xnorm=wavenorm)
        # tailinde = self.coeff['TAILINDE'].eval(ispec, xnorm=wavenorm)
        tailxsca = self.get_coeff('TAILXSCA', ispec, wavelength, iwave)
        tailysca = self.get_coeff('TAILYSCA', ispec, wavelength, iwave)
        tailamp  = self.get_coeff('TAILAMP',  ispec, wavelength, iwave)
        tailcore = self.get_coeff('TAILCORE', ispec, wavelength, iwave)
        tailinde = self.get_coeff('TAILINDE', ispec, wavelength, iwave)

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
        xfunc1 = [self._pgh(xccd, i, x, sigma=sigx1) for i in range(degx1+1)]
        yfunc1 = [self._pgh(yccd, i, y, sigma=sigy1) for i in range(degy1+1)]

        #- Create core PSF image
        core1 = np.zeros((ny, nx))
        spot1 = np.empty_like(core1)
        for i in range(degx1+1):
            for j in range(degy1+1):
                # c1 = self.coeff['GH-{}-{}'.format(i,j)].eval(ispec, xnorm=wavenorm)
                c1 = self.get_coeff('GH-{}-{}'.format(i,j), ispec, wavelength, iwave)
                outer(yfunc1[j], xfunc1[i], out=spot1)
                core1 += c1 * spot1

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
        xc = self.coeff['X'].eval(ispec, wavelength)
        yc = self.coeff['Y'].eval(ispec, wavelength)

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

