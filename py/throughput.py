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
import scipy.linalg
from specter import util

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
        
        #- Flux -> photons conversion constant
        #-         h [erg s]      * c [m/s]      * [1e10 A/m] = [erg A]
        self._hc = 6.62606957e-27 * 2.99792458e8 * 1e10
        
        if fiberinput is not None:
            if isinstance(fiberinput, float):
                self._fiberinput = N.ones(len(wave)) * fiberinput
            else:
                self._fiberinput = N.copy(fiberinput)
        else:
            self._fiberinput = N.ones(len(wave))
    
    @property
    def fiberarea(self):
        return N.pi * self.fiberdia**2 / 4.0
        
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
        
    def _throughput(self, objtype=ObjType.STAR, airmass=1.0):
        """
        Returns system throughput for this object type and airmass
        at the native wavelengths of this Throughput object (self._wave)
        
        objtype may be any of the ObjType enumerated types
            CALIB : atmospheric extinction and fiber input losses not applied
            SKY   : fiber input losses are not applied
            other : all throughput losses are applied
        """
        objtype = objtype.strip().upper()
        
        Tatm = 10**(-0.4*airmass*self._extinction)
        if objtype == ObjType.CALIB:
            T = self._thru
        elif objtype == ObjType.SKY:
            T = self._thru * Tatm
        else:
            T = self._thru * Tatm * self._fiberinput
            
        return T

    def __call__(self, wavelength, objtype=ObjType.STAR, airmass=1.0):
        """
        Returns system throughput at requested wavelength(s)
        
        objtype may be any of the ObjType enumerated types
            CALIB : atmospheric extinction and fiber input losses not applied
            SKY   : fiber input losses are not applied
            other : all throughput losses are applied
        """

        T = self._throughput(objtype=objtype, airmass=airmass)
        return N.interp(wavelength, self._wave, T, left=0.0, right=0.0)

    def thru(self, *args, **kwargs):
        """
        same as calling self(*args, **kwargs)
        """
        return self(*args, **kwargs)

    def photons(self, wavelength, flux, units="erg/s/cm^2/A", \
                objtype="STAR", exptime=None, airmass=1.0):
        """
        Returns photons per bin given input flux vs. wavelength,
        flux units, object type, exposure time, and airmass.
        
        NOTE: Returns raw photons per wavelength bin,
              *not* photons/A sampled at these wavelengths

        Inputs
        ------
        wavelength : input wavelength array in Angstroms
        flux       : input flux; same length as `wavelength`
        units      : units of `flux`
          * Treated as delta functions at each given wavelength:
            - "photons"
            - "erg/s/cm^2"
          * Treated as function values to be multipled by bin width:
            - "photon/A"
            - "erg/s/cm^2/A"
            - "erg/s/cm^2/A/arcsec^2"

        "photon" and "photon/A" are as observed by CCD and are returned
        without applying throughput terms.

        Optional Inputs
        ---------------
        objtype : string, optional; object type for Throughput object.
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
        May 2013
        """
        
        #- Wavelength bin size
        dw = N.gradient(wavelength)
        
        #- Standardize units; allow some sloppiness
        units = units.strip()  #- FITS pads short strings with spaces (!)
        units = units.replace("ergs", "erg")
        units = units.replace("photons", "photon")
        units = units.replace("Angstroms", "A")
        units = units.replace("Angstrom", "A")
        units = units.replace("**", "^")
        units = units.replace("cm2", "cm^2")
        
        #- Check for units prefactor like "1e-17 erg/s/cm^2/A"
        scale = 1.0
        tmp = units.split()
        if len(tmp) == 2:
            try:
                scale = float(tmp[0])
                flux = flux * scale
                units = tmp[1]
            except ValueError:
                raise ValueError, "Non-numeric units scale factor " + tmp[0]
        
        #- Default exposure time
        if exptime is None:
            exptime = self.exptime
        
        #- Input photons; return photons per bin (not photons per Angstrom)
        if units == "photon":
            return flux
        elif units == "photon/A":
            return flux * dw

        #- Sanity check on units
        if not units.startswith('erg'):
            raise ValueError, "Unrecognized units " + units

        #- If we got here, we need to apply throughputs
        flux = self.apply_throughput(wavelength, flux,
                                 objtype=objtype, airmass=airmass)
        
        #- Convert to photons
        phot = flux * wavelength / self._hc
        
        #- erg/s/cm^2/A (i.e. Flambda, astronomical object)
        if units == "erg/s/cm^2/A":
            return phot * exptime * self.effarea * dw
            
        #- erg/s/cm^2/A/arcsec^2 (e.g. sky)
        elif units == "erg/s/cm^2/A/arcsec^2":
            return phot * exptime * self.effarea * dw * self.fiberarea
            
        #- erg/s/cm^2 (not per A; flux delta functions at given wavelengths)
        elif units == "erg/s/cm^2":
            return phot * exptime * self.effarea

        #- erg/s/cm^2/arcsec^2 (not per A; intensity delta functions)
        elif units == "erg/s/cm^2/arcsec^2":
            return phot * exptime * self.effarea * self.fiberarea
            
        else:
            raise ValueError, "Unrecognized units " + units
        
    def apply_throughput(self, wavelength, flux, objtype="STAR", airmass=1.0):
        """
        Returns flux array with throughputs applied for given
        objtype and airmass.
        """

        #- Find function which integrates to each bin flux
        y = util.model_function(wavelength, flux)
        
        #- Apply throughput, sampling at both input wavelengths and
        #- at native throughput wavelengths.
        ww = N.concatenate( (wavelength, self._wave) )
        yflux = N.interp(ww, wavelength, y)
        t = self.thru(ww, objtype=objtype, airmass=airmass)
        yflux *= t
        
        #- HACK: Only apply throughput to positive fluxes
        # ii = N.where(yflux>0)
        # yflux[ii] *= t[ii]
        
        #- Resample back to input wavelength binning
        edges = util.get_bin_edges(wavelength)
        binwidths = N.diff(edges)
        ii = N.argsort(ww)
        binnedflux = util.trapz(edges, ww[ii], yflux[ii]) / binwidths
        return binnedflux

