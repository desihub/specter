from psf import PSF
from spotgrid import SpotGridPSF
from pixpsf import PixPSF
from monospot import MonoSpotPSF
from monospot import MonoSpotPSF
from gausshermite import GaussHermitePSF

def load_psf(filename):
    """
    Load a PSF fits file, using the "PSFTYPE" header keyword to determine
    which specific subclass to use.
    """
    import fitsio
    hdr = fitsio.read_header(filename)
    
    if hdr['PSFTYPE'] == 'SPOTGRID':
        return SpotGridPSF(filename)
    if hdr['PSFTYPE'] == 'MONOSPOT':
        return MonoSpotPSF(filename)
    elif hdr['PSFTYPE'].strip() == 'PCA-PIX':
        return PixPSF(filename)
    elif hdr['PSFTYPE'].strip() == 'GAUSS-HERMITE':
        return GaussHermitePSF(filename)
    else:
        print "Unknown PSFTYPE", hdr['PSFTYPE']
        return PSF(filename)
