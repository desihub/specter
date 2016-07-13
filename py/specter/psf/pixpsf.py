#!/usr/bin/env python

"""
Pixelated 2D PSF

David Schlegel & Stephen Bailey, Summer 2011
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from astropy.io import fits
from specter.psf import PSF
from specter.util import sincshift

#- Turn off complex -> real warnings in sinc interpolation
import warnings 
try:
	warnings.simplefilter("ignore", np.ComplexWarning)
except   AttributeError:
	pass
	
class PixPSF(PSF):
    """
    Pixelated PSF[ny, nx] = Ic[ny,nx] + x*Ix[ny, nx] + y*Iy[ny,nx] + ...
    """
    def __init__(self, filename):
        """
        Loads pixelated PSF parameters from a fits file
        """
        #- Use PSF class to Load Generic PSF values (x, y, wavelength, ...)
        PSF.__init__(self, filename)
        
        #- Additional headers are a custom format for the pixelated psf
        fx = fits.open(filename, memmap=False)
        self.nexp     = fx[3].data.view(np.ndarray)  #- icoeff xexp yexp
        self.xyscale  = fx[4].data.view(np.ndarray)  #- ifiber igroup x0 xscale y0 yscale
        self.psfimage = fx[5].data.view(np.ndarray)  #- [igroup, icoeff, iy, ix]
        fx.close()
                
    def _xypix(self, ispec, wavelength):
        """
        Evaluate PSF for a given spectrum and wavelength
        
        returns xslice, yslice, pixels[yslice, xslice]
        """
        #- Get fiber group and scaling factors for this spectrum
        igroup = self.xyscale['IGROUP'][ispec]
        x0     = self.xyscale['X0'][ispec]
        xscale = self.xyscale['XSCALE'][ispec]
        y0     = self.xyscale['Y0'][ispec]
        yscale = self.xyscale['YSCALE'][ispec]
        
        #- Get x and y centroid
        x, y = self.xy(ispec, wavelength)
        
        #- Rescale units
        xx = xscale * (x - x0)
        yy = yscale * (y - y0)
        
        #- Generate PSF image at (x,y)
        psfimage = np.zeros(self.psfimage.shape[2:4])
        for i in range(self.psfimage.shape[1]):
            nx = self.nexp['XEXP'][i]
            ny = self.nexp['YEXP'][i]
            psfimage += xx**nx * yy**ny * self.psfimage[igroup, i]
                                
        #- Sinc Interpolate
        dx = x - int(round(x))
        dy = y - int(round(y))
        psfimage = sincshift(psfimage, dx, dy)
        
        #- Check boundaries
        ny, nx = psfimage.shape
        ix = int(round(x))
        iy = int(round(y))
        
        xmin = max(0, ix-nx//2)
        xmax = min(self.npix_x, ix+nx//2+1)
        ymin = max(0, iy-ny//2)
        ymax = min(self.npix_y, iy+ny//2+1)
                
        if ix < nx//2:
            psfimage = psfimage[:, nx//2-ix:]
        if iy < ny//2:
            psfimage = psfimage[ny//2-iy:, :]
        
        if ix+nx//2+1 > self.npix_x:
            dx = self.npix_x - (ix+nx//2+1)
            psfimage = psfimage[:, :dx]
            
        if iy+ny//2+1 > self.npix_y:
            dy = self.npix_y - (iy+ny//2+1)
            psfimage = psfimage[:dy, :]
        
        xslice = slice(xmin, xmax)
        yslice = slice(ymin, ymax)

        #- Clip negative values
        psfimage[psfimage<0] = 0.0

        #- Normalize integral to 1.0
        psfimage /= psfimage.sum()
        
        assert xslice.stop-xslice.start == psfimage.shape[1]
        assert yslice.stop-yslice.start == psfimage.shape[0]

        return xslice, yslice, psfimage

