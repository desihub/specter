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
wave1D = dwave*N.arange(nwave)+6000.0
wave2D = N.tile(wave1D, nspec).reshape(nspec, nwave)
waves = [wave1D, wave2D, None]
loglams = [N.log10(wave1D), N.log10(wave2D)]

targetinfo = dict(OBJTYPE = N.array(["STAR", ] * nspec))

#- Generate sequential filenames
n = -1
def get_next_filename():
    global n
    n += 1
    return "data/spec-%03d.fits" % n

#- Wrapper to write a spectrum file as image HDUs
#- Returns name of file written
def write_imgspec(flux, wave, loglam, objtype):
    hdr = dict(FLUXUNIT='erg/s/cm^2/A')
    if wave is None:
        if loglam:
            hdr['LOGLAM'] = loglam
            hdr['CRVAL1'] = N.log10(wave1D[0])
            hdr['CDELT1'] = N.log10(1+dwave/wave1D[0])
        else:
            if loglam is not None:
                hdr['LOGLAM'] = loglam
            hdr['CRVAL1'] = wave1D[0]
            hdr['CDELT1'] = dwave

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

    return filename
    
#- Wrapper to write a spectrum file as a binary table
#- Returns name of file written
def write_tblspec(flux, wave, loglam, objtype):
    data = dict(flux=flux)    
    hdr = dict(FLUXUNIT='erg/s/cm^2/A')
    if wave is None:
        if loglam:
            hdr['LOGLAM'] = loglam
            hdr['CRVAL1'] = N.log10(wave1D[0])
            hdr['CDELT1'] = N.log10(1+dwave/wave1D[0])
        else:
            if loglam is not None:
                hdr['LOGLAM'] = loglam
            hdr['CRVAL1'] = wave1D[0]
            hdr['CDELT1'] = dwave
    else:
        if loglam is None or not loglam:
            data['wavelength'] = wave
        else:
            data['loglam'] = wave

    if type(objtype) == str:
        hdr['OBJTYPE'] = objtype
    else:
        data['objtype'] = objtype['OBJTYPE']
    
    filename = get_next_filename()
    fitsio.write(filename, data, header=hdr, clobber=True)

    return filename

for flux in (flux1D, flux2D):
    for wave in (wave1D, wave2D, None):
        if wave is not None and (wave.ndim > flux.ndim):
            continue
            
        for loglam in (None, False, True):
            if loglam and wave is not None:
                xwave = N.log10(wave)
            else:
                xwave = wave
                
            for objtype in ('STAR', targetinfo):
                if (flux.ndim == 1) and (objtype == targetinfo):
                    continue
                filename = write_imgspec(flux, xwave, loglam, objtype)
                wdim = wave.ndim if wave is not None else 0
                if type(objtype) == dict:
                    objtype_str = ' %-6s' % (objtype['OBJTYPE'][0]+'*',)
                else:
                    objtype_str = ' %-6s' % objtype
                
                print filename, 'image', objtype_str, loglam, wdim, flux.ndim
                if (wave is not None) and (flux.ndim == wave.ndim):
                    filename = write_tblspec(flux, xwave, loglam, objtype)
                    print filename, 'table', objtype_str, loglam, wdim, flux.ndim

