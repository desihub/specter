#!/usr/bin/env python

"""
A class for tracking throughput

Hacks:
  - Doesn't support  spatial variation of input sources
  - Doesn't support per-fiber throughput
  - Do I really want to impose a clumsy ObjType ENUM?
  
How to handle fiber size and sky units?
"""

import sys
import os
import numpy as N
import fitsio

#- ObjType enum
class ObjType:
    STAR   = 'STAR'
    GALAXY = 'GALAXY'
    QSO    = 'QSO'
    SKY    = 'SKY'
    CALIB  = 'CALIB'

def throughput(filename):
    """
    Create Throughput object from FITS file with EXTNAME=THROUGHPUT HDU
    """
    thru = fitsio.read(filename, 'THROUGHPUT', lower=True)
    hdr =  fitsio.read_header(filename, 'THROUGHPUT')
    
    if 'wavelength' in thru.dtype.names:
        w = thru['wavelength']
    elif 'loglam' in thru.dtype.names:
        w = 10**thru['loglam']
    else:
        raise ValueError, 'throughput must include wavelength or loglam'
        
    return Throughput(w, thru['throughput'], thru['extinction'], thru['fiberinput'])

class Throughput:
    def __init__(self, wave, throughput, extinction,
        exptime, effarea, fiberdia, fiberinput=None):
        """
        Create Throughput object
        
        Inputs
        ------
        wave : wavelength array [Angstroms]
        throughput : array of system throughput for elements which apply to
            all types of sources, e.g. mirrors, lenses, CCD efficiency
        extinction : atmospheric extinction array [mags/airmass]
        
        exptime: float, default exposure time [sec]
        effarea: float, primary mirror effective area [cm^2]
        fiberdia: float, fiber diameter [arcsec]
        
        Optional Inputs
        ---------------
        fiberinput : float or array; geometric loss from finite sized
            fiber input.  Default to no loss.
            
            
        Notes
        -----
        fiberinput is a placeholder, since it really depends upon the
        spatial extent of the object and the seeing.
        """
        self._wave = N.copy(wave)
        self._thru = N.copy(throughput)
        self._extinction  = N.copy(extinction)
        
        self.exptime = exptime
        self.effarea = effarea
        self.fiberdia = fiberdia
        
        if fiberinput is not None:
            if isinstance(fiberinput, float):
                self._fiberinput = N.ones(len(wave)) * fiberinput
            else:
                self._fiberinput = N.copy(fiberinput)
        else:
            self._fiberinput = N.ones(len(wave))
        
    def __call__(self, wavelength, objtype=ObjType.STAR, airmass=1.0):
        """
        Returns system throughput at requested wavelength(s)
        
        objtype may be any of the ObjType enumerated types
            CALIB : atmospheric extinction and fiber input losses not applied
            SKY   : fiber input losses are not applied
            other : all throughput losses are applied
        """
        objtype = objtype.upper()
        
        Tatm = 10**(-0.4*airmass*self._extinction)
        if objtype == ObjType.CALIB:
            T = self._thru
        elif objtype == ObjType.SKY:
            T = self._thru * Tatm
        else:
            T = self._thru * Tatm * self._fiberinput
            
        return N.interp(wavelength, self._wave, T)

    def thru(self, *args, **kwargs):
        """
        same as calling self(*args, **kwargs)
        """
        return self(*args, **kwargs)

    def photons(self, wavelength, flux, units, \
                objtype="STAR", exptime=None, airmass=1.0):
        """
        Returns photons observed by CCD given input flux vs. wavelength,
        flux units, and object type, using this throughput model.

        Inputs
        ------
        wavelength : input wavelength array in Angstroms
        flux       : input flux or photons; same length as `wavelength`
        units      : units of `flux`
          * Treated as delta functions at each wavelength:
            - "photons"
            - "ergs/s/cm^2"
          * Treated as function values to be multipled by bin width:
            - "photons/A"
            - "photons/A/arcsec^2"
            - "ergs/s/cm^2/A"
            - "ergs/s/cm^2/A/arcsec^2"

        For the per-Angstrom units options ("/A"), the user is responsible for
        providing a fine enough wavelength sampling that multiplying
        by the bin width is a reasonable approximation of the integration.

        Optional Inputs
        ---------------
        objtype : string, optional; object type for Throuput object.
            If objtype is not 'SKY' or 'CALIB', it will be treated as an
            astronomical source with all throughput terms applied.
        exptime : float, optional; exposure time, default self.exptime
        airmass : float, optional, default 1.0

        Returns
        -------
        array of number of photons observed by CCD at each wavelength,
        i.e. not per-Angstrom.

        Stephen Bailey, LBL
        January 2013
        """

        #- Flux -> photons conversion constants
        #-   h [ergs s] * c [m/s] * [1e10 A/m] = [ergs A]
        h_c = 1.05457168e-27 * 2.99792458e-8  * 1e10

        #- Get default exposure time from Throughput object
        if exptime is None:
            exptime = self.exptime

        #- photons * throughput = easy
        if units == "photons":
            return flux * self.thru(wavelength, objtype=objtype, airmass=airmass)

        #- photons/A
        elif units == "photons/A":
            #- photons/A * width(A) -> photons
            photons = flux * N.gradient(wavelength)

            #- Now multiply by throughput
            return photons * self.thru(wavelength, objtype=objtype, airmass=airmass)

        #- photons/A/arcsec^2 (e.g. calibration lamp)
        elif units == "photons/A/arcsec^2":
            phot = flux * self.thru.fiberdia**2
            return photons(wavelength, phot, units='photons/A', \
                           objtype=objtype, exptime=exptime, airmass=airmass)

        #- ergs/s/cm^2/A (e.g. astronomical object)
        elif units == "ergs/s/cm^2/A":
            phot = flux * wavelength / h_c              #- photons/s/cm^2/A
            phot *= thru.effarea                        #- photons/s/A
            phot *= exptime                             #- photons/A
            return photons(wavelength, phot, units='photons/A', \
                           objtype=objtype, exptime=exptime, airmass=airmass)

        #- ergs/s/cm^2/A/arcsec^2 (e.g. sky)
        elif units == "ergs/s/cm^2/A/arcsec^2":
            f = flux * thru.fiberdia**2 / 4.0
            return photons(ispec, wavelength, f, units="ergs/s/cm2/A", \
                           objtype=objtype, exptime=exptime, airmass=airmass)
    
        