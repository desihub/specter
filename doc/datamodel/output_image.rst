===================
Output Image Format
===================

HDU 0 CCDIMAGE
--------------

Image of spectra projected onto the CCD with the PSF, with optional noise.
The readout noise is always Gaussian, but the photon noise could be
Poisson (default) or Gaussian (if --gaussnoise option was used.)

Header keywords:

- SIMDATA = True
- PREPROC = True
- GAIN    = CCD gain in electrons/ADU
- RDNOISE = CCD amplifier readout noise in electrons
- SIMNOISE = "Gaussian", "Poisson", or "None"

HDU 1 IVAR
----------

If noise was added to the CCD, this contains the pixel inverse variance
if the noise was treated as Gaussian (even if the photon shot noise is
Poisson).

HDU 2 PHOTONS
-------------

*NOTE* : The format of this HDU will likely change to match the format
of input spectra (``flux`` plus ``FLUXUNIT`` keyword instead of ``PHOTONS``).
This has not yet been implemented.

Optional HDU with binary table giving the photon spectra projected onto
the CCD after all throughputs were applied.  Extraction code should
produce this before any calibrations.  There is one row per spectrum,
with columns:

- PHOTONS[nwave]
- WAVELENGTH[nwave]

HDU 3 XYWAVE
------------

Optional HDU with x,y vs. wavelength of spectral traces on CCD.
There is one row per spectrum, with columns:

- X[nwave]
- Y[nwave]
- WAVELENGTH[nwave]
