"""
specter.psf
===========
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from .psf import PSF
from .spotgrid import SpotGridPSF
from .pixpsf import PixPSF
from .monospot import MonoSpotPSF
from .monospot import MonoSpotPSF
from .gausshermite import GaussHermitePSF
from .gausshermite2 import GaussHermite2PSF

def load_psf(filename, psftype=None):
    """
    Load a PSF fits file, using the "PSFTYPE" header keyword to determine
    which specific subclass to use.

    psftype : override fits header keyword PSFTYPE
    """
    from astropy.io import fits
    hdr = fits.getheader(filename)

    if psftype is not None:
        hdr['PSFTYPE'] = psftype

    if hdr['PSFTYPE'] == 'SPOTGRID':
        return SpotGridPSF(filename)
    if hdr['PSFTYPE'] == 'MONOSPOT':
        return MonoSpotPSF(filename)
    elif hdr['PSFTYPE'].strip() == 'PCA-PIX':
        return PixPSF(filename)
    elif hdr['PSFTYPE'].strip() == 'GAUSS-HERMITE':
        return GaussHermitePSF(filename)
    elif hdr['PSFTYPE'].strip() == 'GAUSS-HERMITE2':
        return GaussHermite2PSF(filename)
    else:
        print("Unknown PSFTYPE {}".format(hdr['PSFTYPE']))
        return PSF(filename)
