"""
MonoSpotPSF - ...
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from astropy.io import fits
from specter.psf import PSF
from specter.util import LinearInterp2D, rebin_image, sincshift

class MonoSpotPSF(PSF):

    def __init__(self, filename, spot=None, scale=1.0):
        """
        Initialize MonoSpotPSF from input file with optional override of
        which spot[y,x] to use.  If overriding spot, scale gives ratio of
        spot pixel-size to CCD pixel-size.  Must be >1 and evenly divisible
        by spot dimensions.
        
        See specter.psf.PSF for futher details
        """
        #- Use PSF class to Load Generic PSF values (x, y, wavelength, ...)
        PSF.__init__(self, filename)
        
        if spot is None:
            self._spot = fits.getdata(filename, 'SPOT')
            self._scale = fits.getheader(filename, 'SPOT')['SCALE']
        else:
            self._spot = spot.copy()
            self._scale = scale
        
    def _xypix(self, ispec, wavelength):
        """
        Return xslice, yslice, pix for PSF at spectrum ispec, wavelength
        """
        assert 0 <= ispec < self.nspec
        
        xc, yc = self.xy(ispec, wavelength)
        scale = self._scale  #- shorthand

        #- Calculate offset into CCD pixel
        xoffset = int(xc * scale) % scale
        yoffset = int(yc * scale) % scale

        #- Place high res spot into grid aligned with CCD pixels
        ny, nx = self._spot.shape
        A = np.zeros(shape=(ny+scale, nx+scale))
        A[yoffset:yoffset+ny, xoffset:xoffset+nx] = self._spot
        ccdpix = rebin_image(A, scale)

        #- Fractional high-res pixel offset
        #- This can be slow; is it really necessary?
        dxx = ((xc * scale) % scale - xoffset) / scale
        dyy = ((yc * scale) % scale - yoffset) / scale
        ccdpix = sincshift(ccdpix, dxx, dyy)

        #- sinc shift can cause negative ringing, so clip and re-normalize
        ccdpix = ccdpix.clip(0)
        ccdpix /= np.sum(ccdpix)

        #- Find where the [0,0] pixel goes on the CCD 
        xccd = int(xc - ccdpix.shape[1]//2 + 1)
        yccd = int(yc - ccdpix.shape[0]//2 + 1)

        xx = slice(xccd, xccd+ccdpix.shape[1])
        yy = slice(yccd, yccd+ccdpix.shape[0])

        return xx, yy, ccdpix
    

#- Incomplete code for creating without a file
# def __init__(self, x, y, w, spot, scale=1):
#     """
#     Create new MonoSpotPSF
#     
#     Inputs:
#         x[nspec, *] : 2D array of x-centroids on CCD
#         y[nspec, *] : 2D array of y-centroids on CCD
#         w[nspec, *] : 2D array of wavelengths on CCD
#         
#         spot[ny,nx] : 2D spot image to use everywhere on CCD
#         scale : ratio of spot pixel scale to CCD pixel scale.
#                 Must be evenly divisible by ny,nx.
# 
#     Notes:
#     The second dimension of x,y,w is arbitrary, but the resolution
#     should be fine enough that liner interpolation is valid.
#     e.g. numpy.interp(yt, y[i], w[i]) should give the wavelength of
#     spectrum i at CCD row = yt.
#     
#     If spot.shape == (120,120) and scale=10, then 10 pixels of spot
#     correspond to 1 CCD pixel and the spot covers an area of (12x12)
#     CCD pixels.
#     """
# 
#     assert x.shape == y.shape == w.shape
#     assert spot.shape[0] % scale == 0
#     assert spot.shape[1] % scale == 0
#     
#     self._nspec = x.shape[0]
#     self._x = x.copy()
#     self._y = y.copy()
#     self._w = w.copy()
#     self._spot = spot.copy()
#     self._scale = scale
