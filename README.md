specter : 2D CCD pixel-level simulation of multi-object spectrographs

Directories included here:

doc : documentation (data model, overview, ...)
py  : python classes and modules
    py/psf : modeling the spectrograph point spread functions
    py/extract : extract spectra from CCD images given a PSF (TODO)
bin : executable scripts
dev : code under initial development, before it moves into py/
lib : add this to PYTHONPATH to get the live development version

Run bin/specter -h to see command line options.

To use the python toolkit:

    #- Load PSF file for some spectrograph
    from specter.psf import load_psf
    psf = load_psf('psf.fits')
    
    #- Get PSF spot for fiber 10 and wavelength 5000 Angstroms:
    spot = psf.pix(10, 5000)
    
    #- Project flux[nspec, nw] at wavelength[nw] onto the CCD:
    image = psf.project(flux, wavelength)

Stephen Bailey
Fall 2012

