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

def load_throughput(filename):
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
        
    return Throughput(
        wave = w,
        throughput = thru['throughput'],
        extinction = thru['extinction'],
        fiberinput = thru['fiberinput'],
        exptime  = hdr['EXPTIME'],
        effarea  = hdr['EFFAREA'],
        fiberdia = hdr['FIBERDIA'],
        )

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
        fiberinput : float or array; geometric throughput from finite sized
            fiber input.  Default to no loss = 1.0.
            
            
        Notes
        -----
        fiberinput is a placeholder, since it really depends upon the
        spatial extent of the object and the seeing.
        """
        self._wave = N.copy(wave)
        self._thru = N.copy(throughput)
        self._extinction  = N.copy(extinction)
        
        self.exptime = float(exptime)
        self.effarea = float(effarea)
        self.fiberdia = float(fiberdia)
        
        if fiberinput is not None:
            if isinstance(fiberinput, float):
                self._fiberinput = N.ones(len(wave)) * fiberinput
            else:
                self._fiberinput = N.copy(fiberinput)
        else:
            self._fiberinput = N.ones(len(wave))
        
    def extinction(self, wavelength):
        """
        Return atmospheric extinction [magnitudes/airmass]
        evaluated at wavelength (float or array)
        """
        return N.interp(wavelength, self._wave, self._extinction)
    
    def atmospheric_throughput(self, wavelength, airmass=1.0):
        """
        Return atmospheric throughput [0 - 1] at given airmass,
        evaluated at wavelength (float or array)
        """
        ext = self.extinction(wavelength)
        return 10**(-0.4 * airmass * ext)
        
    def fiberinput_throughput(self, wavelength):
        """
        Return fiber input geometric throughput [0 - 1]
        evaluated at wavelength (float or array)
        """
        return N.interp(wavelength, self._wave, self._fiberinput)
        
    def hardware_throughput(self, wavelength):
        """
        Return hardware throughput (optics, fiber run, CCD, but
        not including atmosphere or geometric fiber input losses)
        evaluated at wavelength (float or array)        
        """
        return N.interp(wavelength, self._wave, self._thru, left=0.0, right=0.0)
        
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
            
        return N.interp(wavelength, self._wave, T, left=0.0, right=0.0)

    def thru(self, *args, **kwargs):
        """
        same as calling self(*args, **kwargs)
        """
        return self(*args, **kwargs)

    def photons(self, wavelength, flux, units="erg/s/cm^2/A", \
                objtype="STAR", exptime=None, airmass=1.0):
        """
        Returns photons observed by CCD given input flux vs. wavelength,
        flux units, and object type, using this throughput model.

        Inputs
        ------
        wavelength : input wavelength array in Angstroms
        flux       : input flux or photons; same length as `wavelength`
        units      : units of `flux`
          * Treated as delta functions at each given wavelength:
            - "photons"
            - "erg/s/cm^2"
          * Treated as function values to be multipled by bin width:
            - "photon/A"
            - "erg/s/cm^2/A"
            - "erg/s/cm^2/A/arcsec^2"

        For the per-Angstrom units options ("/A"), the user is responsible for
        providing a fine enough wavelength sampling that multiplying
        by the bin width is a reasonable approximation of the integration.

        "photon" and "photon/A" are as observed by CCD and are returned
        without applying throughput terms.

        Optional Inputs
        ---------------
        objtype : string, optional; object type for Throuput object.
            SKY - atmospheric extinction and telescope+instrument throughput
                    applied, but not fiber input geometric losses.
            CALIB - telescope+instrument throughtput applied, but not
                    atmospheric extinction or fiber input geometric losses.
            Anything else (default) - treated as astronomical object
                    with all throughput terms applied.
        exptime : float, optional; exposure time, default self.exptime
        airmass : float, optional, default 1.0

        Returns
        -------
        array of number of photons observed by CCD at each wavelength,
        i.e. not per-Angstrom.

        Stephen Bailey, LBL
        January 2013
        """

        #- Allow some sloppiness in units
        units = units.strip()  #- FITS pads short strings with spaces (!)
        units = units.replace("ergs", "erg")
        units = units.replace("photons", "photon")
        units = units.replace("Angstroms", "A")
        units = units.replace("Angstrom", "A")
        units = units.replace("**", "^")

        #- check if units have numerical scaling prefix
        #- e.g. "1e-16 erg/s/cm^2/A"
        try:
            pre, xunit = units.split()
            flux = float(pre) * flux
            units = xunit
        except:
            pass

        #- Flux -> photons conversion constants
        #-   h [erg s] * c [m/s] * [1e10 A/m] = [erg A]
        h_c = 6.62606957e-27 * 2.99792458e8 * 1e10

        #- Default exposure time
        if exptime is None:
            exptime = self.exptime

        #- photons (delta functions at given wavelengths)
        if units == "photon":
            return flux

        #- photon/A
        elif units == "photon/A":
            #- photon/A * width(A) -> photons
            phot = flux * N.gradient(wavelength)
            return phot

        #- photon/A/arcsec^2
        elif units == "photon/A/arcsec^2":
            #- photon/A/arcsec^2 * width(A) -> photons/arcsec^2
            phot = flux * N.gradient(wavelength)
            #- multiply by area to get final photons
            phot *= self.effarea
            return phot

        #- erg/s/cm^2 (flux delta functions at given wavelengths)
        elif units == "erg/s/cm^2":
            phot = flux * wavelength / h_c              #- photon/s/cm^2
            phot *= self.effarea                        #- photon/s
            phot *= exptime                             #- photon
            return phot * self.thru(wavelength, objtype=objtype, airmass=airmass)
            
        elif units == "erg/s/cm^2/arcsec^2":
            phot = flux * wavelength / h_c              #- photon/s/cm^2/arcsec^2
            phot *= self.effarea                        #- photon/s/arcsec^2
            phot *= exptime                             #- photon/arcsec^2
            phot *= N.pi * self.fiberdia**2 / 4.0       #- photon
            return phot * self.thru(wavelength, objtype=objtype, airmass=airmass)
            
        #- erg/s/cm^2/A (e.g. astronomical object)
        elif units == "erg/s/cm^2/A":
            phot = flux * wavelength / h_c              #- photon/s/cm^2/A
            phot *= self.effarea                        #- photon/s/A
            phot *= exptime                             #- photon/A
            phot *= N.gradient(wavelength)              #- photons
            return phot * self.thru(wavelength, objtype=objtype, airmass=airmass)

        #- erg/s/cm^2/A/arcsec^2 (e.g. sky)
        elif units == "erg/s/cm^2/A/arcsec^2":
            f = flux * N.pi * self.fiberdia**2 / 4.0    #- erg/s/cm^2/A
            return self.photons(wavelength, f, units="erg/s/cm^2/A", \
                           objtype=objtype, exptime=exptime, airmass=airmass)

        else:
            raise ValueError, "Unrecognized units " + units
    
        