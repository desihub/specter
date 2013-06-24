"""
Routines to standarize specter I/O.  This may be refactored into a
separate directory as this grows.

Stephen Bailey, LBL
January 2013
"""

import fitsio
import numpy as N

def read_simspec(filename):
    """
    Read an input simulated spectrum file, parse the various format
    options, and return a standardized spectrum dictionary for use.
    
    Returns a dictionary with keys flux, wavelength, units, objtype
    """
    
    fx = fitsio.FITS(filename)
    is_image = fx[0].has_data()
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
    fx = fitsio.FITS(filename)
    flux = fx[0].read()
    header = fx[0].read_header()
    if 'wavelength' in fx.hdu_map:
        w = fx['WAVELENGTH'].read()
    elif 'loglam' in fx.hdu_map:
        w = 10**fx['loglam'].read()
    else:
        nwave = flux.shape[-1]
        w = header['CRVAL1'] + N.arange(nwave) * header['CDELT1']
        if 'LOGLAM' in header and header['LOGLAM']:
            w = 10**w
        #- DC-FLAG is deprecated, but still support it
        elif 'DC-FLAG' in header and header['DC-FLAG']:
            w = 10**w

    #- Convert wavelength to 2D if needed
    if flux.ndim == 2 and w.ndim == 1:
        nspec, nwave = flux.shape
        w = N.tile(w, nspec).reshape(nspec, nwave)
        
    #- Get object type or use default
    if 'OBJTYPE' in header:
        objtype = header['OBJTYPE']
    else:
        objtype = 'STAR'
        
    #- Get object type (CALIB, SKY, etc.)
    return dict(flux=flux, wavelength=w,
                units=header['FLUXUNIT'], objtype=objtype)

def read_simspec_table(filename):
    """
    Read an input simulated spectrum file formatted as a FITS binary table.
    
    Returns a dictionary with keys flux, wavelength, units, objtype
    """
    spectra = fitsio.read(filename, 1, lower=True).view(N.recarray)
    header = fitsio.read_header(filename, 1)
        
    #- Extract wavelength in various formats
    if 'wavelength' in spectra.dtype.names:
        w = spectra.wavelength
    elif 'wave' in spectra.dtype.names:
        w = spectra.wave
    elif 'loglam' in spectra.dtype.names:
        w = 10**spectra.loglam
    else:
        nwave = spectra.flux.shape[-1]
        w = header['CRVAL1'] + N.arange(nwave) * header['CDELT1']
        if 'LOGLAM' in header and header['LOGLAM']:
            w = 10**w
        #- DC-FLAG is deprecated, but still support it
        elif 'DC-FLAG' in header and header['DC-FLAG']:
            w = 10**w

    #- Convert wavelength to 2D if needed
    if spectra.flux.ndim == 2 and w.ndim == 1:
        nspec, nwave = spectra.flux.shape
        w = N.tile(w, nspec).reshape(nspec, nwave)
        
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
    else:
        print >> sys.stderr, 'WARNING: using default flux units of erg/s/cm^2/A'
        units = 'erg/s/cm^2/A'
        
    #- return results
    return dict(flux=spectra.flux, wavelength=w, units=units, objtype=objtype)
    