"""
specter
=======

A toolkit for simulating multi-object spectrographs.

The specter toolkit is centered around a generic PSF object, which
describes how a given fiber (spectrum) at a particular wavelength
is projected onto the detector CCD.  The base class :class:`specter.psf.psf.PSF`
defines the interface for any PSF object.  Specific instruments
(*e.g.* BOSS or DESI) implement a subclass which provides::

    xslice, yslice, pixels[ny,nx] = PSF.xypix(ispec, wavelength)

to describe how a given spectrum ispec and wavelength project onto the CCD.
"""
import specter.io

from ._version import __version__
