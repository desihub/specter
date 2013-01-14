def photons(wavelength, flux, units, thru, objtype="STAR", exptime=None):
    """
    Returns photons observed by CCD given input flux vs. wavelength,
    flux units, object type, and throughput model.
    
    Inputs
    ------
    wavelength : input wavelength array in Angstroms
    flux       : input flux or photons; same length as `wavelength`
    units      : units of `flux`
      * Treated as delta functions at each wavelength:
        - "photons"
        - "ergs/s/cm^2"
      * Treated as function values to be multipled by bin width:
        - "photons/A"
        - "photons/A/arcsec^2"
        - "ergs/s/cm^2/A"
        - "ergs/s/cm^2/A/arcsec^2"
    thru : Throughput object

    For the per-Angstrom units options ("/A"), the user is responsible for
    providing a fine enough wavelength sampling that multiplying
    by the bin width is a reasonable approximation of the integration.

    Optional Inputs
    ---------------
    objtype : string, optional; object type for Throuput object.
        If objtype is not 'SKY' or 'CALIB', it will be treated as an
        astronomical source with all throughput terms applied.
    exptime : float, optional; exposure time, default thru.exptime
    
    Returns
    -------
    array of number of photons observed by CCD at each wavelength,
    i.e. not per-Angstrom.
    
    Stephen Bailey, LBL
    January 2013
    """

    #- Flux -> photons conversion constants
    #-   h [ergs s] * c [m/s] * [1e10 A/m] = [ergs A]
    h_c = 1.05457168e-27 * 2.99792458e-8  * 1e10

    #- Get default exposure time from Throughput object
    if exptime is None:
        exptime = thru.exptime

    #- photons * throughput = easy
    if units == "photons":
        return flux * thru.thru(wavelength, objtype=objtype)

    #- photons/A
    elif units == "photons/A":
        #- photons/A * width(A) -> photons
        photons = flux * N.gradient(wavelength)
        
        #- Now multiply by throughput
        return photons * thru.throughput(wavelength, objtype=objtype)

    #- photons/A/arcsec^2 (e.g. calibration lamp)
    elif units == "photons/A/arcsec^2":
        phot = flux * thru.fiberdia**2
        return photons(wavelength, phot, units='photons/A', \
                       objtype=objtype, thru=thru, exptime=exptime)
        
    #- ergs/s/cm^2/A (e.g. astronomical object)
    elif units == "ergs/s/cm^2/A":
        phot = flux * wavelength / h_c              #- photons/s/cm^2/A
        phot *= thru.effarea                        #- photons/s/A
        phot *= exptime                             #- photons/A
        return photons(wavelength, phot, units='photons/A', objtype=objtype, thru=thru, exptime=exptime)
        
    #- ergs/s/cm^2/A/arcsec^2 (e.g. sky)
    elif units == "ergs/s/cm^2/A/arcsec^2":
        f = flux * thru.fiberdia**2 / 4.0
        return photons(ispec, wavelength, f, units="ergs/s/cm2/A", objtype=objtype, thru=thru, psf=psf)
