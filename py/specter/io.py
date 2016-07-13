"""
specter.io
==========

Routines to standarize specter I/O.  This may be refactored into a
separate directory as this grows.

Stephen Bailey, LBL
January 2013
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import os.path
from astropy.io import fits
import numpy as np

def read_image(filename):
    """
    Read (image, ivar, header) from input filename
    """
    fx = fits.open(filename)
    if 'IMAGE' in fx:
        image = fx['IMAGE'].data
        hdr = fx['IMAGE'].header
    else:
        image = fx[0].data
        hdr = fx[0].header

    if 'IVAR' in fx:
        ivar = fx['IVAR'].data
    else:
        ivar = fx[0].data

    fx.close()
    return image, ivar, hdr

def write_spectra(outfile, wave, flux, ivar, resolution, header):
    """
    Write spectra to outfile

    Args:
        wave : 1D[nwave] array of wavelengths
        flux : 2D[nspec, nwave] flux
        ivar : 2D[nspec, nwave] inverse variance of flux
        resolution : 3D[nspec, ndiag, nwave] sparse resolution matrix data
        header : fits header to include in output

    Writes data to outfile in 4 HDUs with EXTNAME FLUX, IVAR, WAVE, RESOLUTION
    """
    assert wave.ndim == 1
    assert flux.ndim == 2
    assert resolution.ndim == 3
    assert wave.shape[0] == flux.shape[1]
    assert flux.shape == ivar.shape
    assert flux.shape[0] == resolution.shape[0]
    assert flux.shape[1] == resolution.shape[2]

    outdir = os.path.dirname(outfile)
    if (outdir != '') and (not os.path.exists(outdir)):
        os.makedirs(outdir)

    hx = fits.HDUList()
    header['EXTNAME'] = 'FLUX'
    hx.append(fits.PrimaryHDU(flux.astype(np.float32), header=header))
    hx.append(fits.ImageHDU(ivar.astype(np.float32), name='IVAR'))
    hx.append(fits.ImageHDU(wave.astype(np.float32), name='WAVELENGTH'))
    hx.append(fits.ImageHDU(resolution.astype(np.float32), name='RESOLUTION'))
    hx.writeto(outfile, clobber=True)


def read_simspec(filename):
    """
    Read an input simulated spectrum file, parse the various format
    options, and return a standardized spectrum dictionary for use.

    Returns a dictionary with keys flux, wavelength, units, objtype
    """

    fx = fits.open(filename)
    is_image = fx[0].data is not None
    fx.close()

    if is_image:
        return read_simspec_image(filename)
    else:
        return read_simspec_table(filename)

def read_simspec_image(filename):
    """
    Read an input simulated spectrum file formatted in multi-HDU FITS images.

    Returns a dictionary with keys flux, wavelength, units, objtype
    """
    fx = fits.open(filename)
    flux = fx[0].data
    header = fx[0].header
    if 'wavelength' in fx:
        w = fx['WAVELENGTH'].data
    elif 'loglam' in fx:
        w = 10**fx['loglam'].data
    else:
        nwave = flux.shape[-1]
        w = header['CRVAL1'] + np.arange(nwave) * header['CDELT1']
        if 'LOGLAM' in header and header['LOGLAM']:
            w = 10**w
        #- DC-FLAG is deprecated, but still support it
        elif 'DC-FLAG' in header and header['DC-FLAG']:
            w = 10**w

    #- Convert wavelength to 2D if needed
    if flux.ndim == 2 and w.ndim == 1:
        nspec, nwave = flux.shape
        w = np.tile(w, nspec).reshape(nspec, nwave)

    #- Get object type or use default
    if 'OBJTYPE' in header:
        objtype = header['OBJTYPE']
    else:
        objtype = 'STAR'

    if 'BUNIT' in header:
        fluxunits = header['BUNIT']
    elif 'FLUXUNIT' in header:
        fluxunits = header['FLUXUNIT']
    else:
        raise ValueError("Unable to determine flux units; need either BUNIT or FLUXUNIT keyword")

    fx.close()
    #- Get object type (CALIB, SKY, etc.)
    return dict(flux=flux, wavelength=w,
                units=fluxunits, objtype=objtype)

def read_simspec_table(filename):
    """
    Read an input simulated spectrum file formatted as a FITS binary table.

    Returns a dictionary with keys flux, wavelength, units, objtype
    """
    spectra, header = fits.getdata(filename, 1, header=True)

    #- Extract wavelength in various formats
    if 'wavelength' in spectra.dtype.names:
        w = spectra.wavelength
    elif 'wave' in spectra.dtype.names:
        w = spectra.wave
    elif 'loglam' in spectra.dtype.names:
        w = 10**spectra.loglam
    else:
        nwave = spectra.flux.shape[-1]
        w = header['CRVAL1'] + np.arange(nwave) * header['CDELT1']
        if 'LOGLAM' in header and header['LOGLAM']:
            w = 10**w
        #- DC-FLAG is deprecated, but still support it
        elif 'DC-FLAG' in header and header['DC-FLAG']:
            w = 10**w

    #- Convert wavelength to 2D if needed
    if spectra.flux.ndim == 2 and w.ndim == 1:
        nspec, nwave = spectra.flux.shape
        w = np.tile(w, nspec).reshape(nspec, nwave)

    #- extract objtype
    if 'objtype' in spectra.dtype.names:
        objtype = spectra.objtype
    else:
        objtype = header['OBJTYPE']

    #- Determine flux units from TUNITnn or FLUXUNIT keywords
    key = 'TUNIT%d' % (spectra.dtype.names.index('flux')+1, )
    if key in header:
        units = header[key].strip()
    elif 'FLUXUNIT' in header:
        units = header['FLUXUNIT'].strip()
    elif 'BUNIT' in header:
        units = header['BUNIT'].strip()
    else:
        print('WARNING: using default flux units of erg/s/cm^2/A', file=sys.stderr)
        units = 'erg/s/cm^2/A'

    #- return results
    return dict(flux=spectra.flux, wavelength=w, units=units, objtype=objtype)
