#!/usr/bin/env python

"""
Write a suite of input spectra files with different formats, units,
1D vs. 2D, etc. to use with testing.
"""

import sys
import os
import numpy as N
import fitsio

nspec = 10
nwave = 20

#- Test data to write
flux1D = N.random.normal(size=nwave)
flux2D = N.random.normal(size=(nspec, nwave))
fluxes = [flux1D, flux2D]

dwave = 1.0
wave1D = dwave*N.arange(nwave)+1.0
wave2D = N.tile(wave1D, nspec).reshape(nspec, nwave)
waves = [wave1D, wave2D, None]
loglams = [N.log10(wave1D), N.log10(wave2D)]
wavetypes = ['WAVELENGTH', 'LOGLAM', 'KEYWORD']

targetinfo = dict(OBJTYPE = N.array(["STAR", ] * nspec))

#- Generate sequential filenames
n = -1
def get_next_filename():
    global n
    n += 1
    return "data/spec-%03d.fits" % n

#- Wrapper to write a spectrum file as image HDUs
def write_imgspec(flux, wave, loglam, objtype):
    hdr = dict(FLUXUNIT='erg/s/cm^2/A')
    if wave is None:
        hdr['CRVAL1'] = wave1D[0]
        hdr['CDELT1'] = dwave
        if loglam is not None:
            hdr['LOGLAM'] = loglam

    if type(objtype) == str:
        hdr['OBJTYPE'] = objtype
    
    filename = get_next_filename()
    fitsio.write(filename, flux, header=hdr, clobber=True)
    if wave is not None:
        if loglam:
            fitsio.write(filename, wave, extname='LOGLAM')
        else:
            fitsio.write(filename, wave, extname='WAVELENGTH')
            
    if type(objtype) != str:
        fitsio.write(filename, objtype, extname='TARGETINFO')
    
#- Wrapper to write a spectrum file as a binary table
def write_tblspec(flux, wave, loglam, objtype):
    data = dict(flux=flux)    
    hdr = dict(FLUXUNIT='erg/s/cm^2/A')
    if wave is None:
        hdr['CRVAL1'] = wave1D[0]
        hdr['CDELT1'] = dwave
        if loglam is not None:
            hdr['LOGLAM'] = loglam
    else:
        if loglam is None or loglam:
            data['wavelength'] = wave
        else:
            data['loglam'] = wave

    if type(objtype) == str:
        hdr['OBJTYPE'] = objtype
    else:
        data['objtype'] = objtype['OBJTYPE']
    
    filename = get_next_filename()
    fitsio.write(filename, data, header=hdr, clobber=True)
    
    
for flux in (flux1D, flux2D):
    for wave in (wave1D, wave2D, None):
        if wave is not None and (wave.ndim > flux.ndim):
            continue
            
        for loglam in (None, True, False):
            for objtype in ('STAR', targetinfo):
                if (flux.ndim == 1) and (objtype == targetinfo):
                    continue
                write_imgspec(flux, wave, loglam, objtype)
                if (wave is not None) and (flux.ndim == wave.ndim):
                    write_tblspec(flux, wave, loglam, objtype)

