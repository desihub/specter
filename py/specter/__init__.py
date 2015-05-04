"""
A toolkit for simulating multi-object spectrographs

The specter toolkit is centered around a generic PSF object, which
describes how a given fiber (spectrum) at a particular wavelength
is projected onto the detector CCD.  The base class specter.psf.PSF
defines the interface for any PSF object.  Specific instruments
(e.g. BOSS, or DESI) implement a subclass which provides:

xslice, yslice, pixels[ny,nx] = PSF.xypix(ispec, wavelength)

to describe how a given spectrum ispec and wavelength project onto the CCD.
"""
import specter.io

__version__ = "0.1.dev1"

def version():
    """Return a version including a git revision string"""
    import os
    from os.path import dirname
    import subprocess
    
    gitdir = dirname(dirname(dirname(__file__)))+'/.git'
    os.environ['GIT_DIR'] = gitdir
    cmd = 'git rev-parse --short=6 HEAD'
    try:
        ver = subprocess.check_output(cmd.split())
        return __version__+'-'+ver.strip()
    except:
        return __version__
        