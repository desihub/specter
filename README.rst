=======
Specter
=======

|Actions Status| |Coveralls Status| |Documentation Status|

.. |Actions Status| image:: https://github.com/desihub/specter/workflows/CI/badge.svg
    :target: https://github.com/desihub/specter/actions
    :alt: GitHub Actions CI Status

.. |Coveralls Status| image:: https://coveralls.io/repos/github/desihub/specter/badge.svg
    :target: https://coveralls.io/github/desihub/specter
    :alt: Test Coverage Status

.. |Documentation Status| image:: https://readthedocs.org/projects/desi-specter/badge/?version=latest
    :target: https://desi-specter.readthedocs.io/en/latest/
    :alt: Documentation Status

Introduction
------------

A 2D CCD pixel-level simulation and extraction toolkit for multi-object spectrographs.

Directories included here:

doc/
    Documentation (data model, overview, ...)

py/specter/
    Python classes and modules

    py/specter/psf/
        modeling the spectrograph point spread functions

    py/specter/extract/
        extract spectra from CCD images given a PSF

bin/
    executable scripts

dev/
    Code under initial development, before it moves into py/


Basic Concepts
--------------

Specter is designed to simulate astronomical multi-object spectrographs,
and extract spectra from the 2D CCD images.  The user provides a model of
the point-spread-function (PSF) of how light from spectrum ``i`` at wavelength
``w`` is dispersed on the CCD::

    pixels, xx, yy = psf.pix(i, w)

where xx and yy are the slices defining the CCD pixel coordinates where
those pixels appear on the CCD image.  doc/datamodel/psf.md lists multiple
file formats for parameterizing different PSFs, or users can subclass
specter.psf.PSF to provide their own parameterization.

Using this PSF model, specter provides pixel-level simulations of the CCD
images produced by the spectrograph given input spectra.

The inverse of simulating photons on the CCD is extracting the photons back
off of the CCD to provide a model of spectra given an input image and noise
model of that image.  Specter uses the PSF model and the
"spectroperfectionism" methodology of Bolton & Schlegel 2009.

Command Line Tools
------------------

Run ``bin/specter -h`` to see command line options for the CCD pixel-level
simulator.  It accepts input spectra, a spectrograph point-spread-function
(PSF) model, and generates output CCD images.

``bin/exspec`` performs 2D PSF extractions given an input image, image noise
model, and PSF.  It subdivides the problem into overlapping regions and
then reassembles them, thus making the O(N^2) problem tractable by doing
many small extractions instead of a single large one.

Python Tools
------------

Run ``python setup.py build_sphinx`` to see python library documentation for both the
simulator and extractor.

Example usage of the Python toolkit::

    #- Load PSF file for some spectrograph
    from specter.psf import load_psf
    psf = load_psf('psf.fits')

    #- Get PSF spot for fiber 10 and wavelength 5000 Angstroms:
    spot = psf.pix(10, 5000)

    #- Project photons[nspec, nw] at wavelength[nw] onto the CCD:
    image = psf.project(photons, wavelength)

    #- Extract spectra 0:10
    from specter.extract import ex2d
    xphot, ivar, R = ex2d(image, image_ivar, psf, specrange=[0,10], \
                          wavelengths=wavelength)

Full Documentation
------------------

Please visit `specter on Read the Docs`_

.. _`specter on Read the Docs`: https://desi-specter.readthedocs.io/en/latest/

License
-------

specter is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.
