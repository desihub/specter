#!/usr/bin/env python

"""
IN PROGRESS; Fit for redshift given a template
"""

import sys
import os
import numpy as N
from scipy.ndimage.filters import gaussian_filter
from specter.util import resample

def zchi2(spectra, wt, template, z):
    wz = wt*(1+z)
    chi2 = 0.0
    for (wave, flux, ivar, psf, ispec) in spectra:
        #- Wavelength at each CCD row
        wpsf = psf.wavelength(ispec, y=N.arange(psf.npix_y))

        #- Template pixel size in Angstroms at those wavelengths
        dw = N.interp(wpsf, wz, N.gradient(wz))        
        
        #- wavelength dispersion in units of template pixels
        wdisp = N.median( psf.wdisp(wpsf) / dw )
        
        #- Convolve with resolution and then sample with input spectrum binning
        model = resample(wz, gaussian_filter(template, wdisp_pix), xnew=wave)

        chi2 += N.sum( (flux-model)**2 * ivar )
        
    return chi2

def zfit(spectra, wt, template, zmin=0.0, zmax=2.0, dlogz=0.001):
    """
    spectra: array of spectra, each of which is a (wave, flux, ivar, psf, ispec) tuple
    wt : wavelength array of template spectrum
    template : template flux
    """
    
    zz = 10**N.arange(N.log10(1.0+zmin), N.log10(1.0+zmax), dlogz)
    chi2 = N.zeros(len(zz))
    for i in range(zz):
        print z[i], i, len(zz)
        chi2[i] = zchi2(spectra, wt, template, zz[i])
        
    return zz, chi2
    
    
    
    
    