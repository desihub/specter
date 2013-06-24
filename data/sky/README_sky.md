Sky Data Files
--------------

This directory contains a few useful files for building up input models
and throughputs.

### ZenExtinct-KPNO.fits ###

HDU 1 has a binary table with columns `WAVELENGTH` [A] and `EXTINCTION`
[mag/airmass] for KPNO.  Below 5000 A the curve is from kpnoextinct.dat;
above 5000 A is from Arjun Dey.  The data were merged by Nick
Mostek for bbspecsim.

### sky-uves.fits ###

HDU 1 has a binary table with columns `wavelength` [A] and `flux`
[1e-17 erg/s/cm^2/A].  This is a sky spectrum from 3400-10600 A stitched
together from the individual spectra described at:

http://www.eso.org/observing/dfo/quality/UVES/pipeline/sky_spectrum.html

