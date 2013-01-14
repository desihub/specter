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

#- ObjType enum
class ObjType:
    STAR   = 'STAR'
    GALAXY = 'GALAXY'
    QSO    = 'QSO'
    SKY    = 'SKY'
    CALIB  = 'CALIB'
        
class Throughput:
    def __init__(self, wave, throughput, extinction, fiberinput):
        """
        Create Throughput object
        
        All inputs are 1D arrays of the same length:
            wave : wavelength array in Angstroms
            throughput : system throughput of elements which apply to
                all types of sources, e.g. mirrors, lenses, CCD efficiency
            extinction : atmospheric extinction in mags per airmass
            fiberinput : geometric loss from finite sized fiber input.
            
        fiberinput is a placeholder, since it really depends upon the
        spatial extent of the object and the seeing.
        """
        self._wave = N.copy(wave)
        self._thru = N.copy(throughput)
        self._extinction  = N.copy(extinction)
        self._fiberinput = N.copy(fiberinput)
        
    def throughput(self, wavelength, objtype=ObjType.STAR, airmass=1.0):
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
        
    def photons(self, wavelength, flux, objtype=ObjType.STAR, airmass=1.0):
        """
        """
        raise NotImplementedError
        T = self.throughput(wavelength, objtype=objtype, airmass=airmass)
        pass  ### TODO ###
        