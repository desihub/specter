from psf import PSF
from spotgrid import SpotGridPSF

def load_psf(filename):
    import fitsio
    hdr = fitsio.read_header(filename)
    
    if hdr['PSFTYPE'] == 'SPOTGRID':
        return SpotGridPSF(filename)
    else:
        print "Unknown PSFTYPE", hdr['PSFTYPE']
        return PSF(filename)
