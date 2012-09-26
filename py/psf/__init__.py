from psf import PSF
from spotgrid import SpotGridPSF

def load_psf(filename):
    """
    Load a PSF fits file, using the "PSFTYPE" header keyword to determine
    which specific subclass to use.
    """
    import fitsio
    hdr = fitsio.read_header(filename)
    
    if hdr['PSFTYPE'] == 'SPOTGRID':
        return SpotGridPSF(filename)
    else:
        print "Unknown PSFTYPE", hdr['PSFTYPE']
        return PSF(filename)
