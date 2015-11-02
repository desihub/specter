=================
Throughput Format
=================

HDU EXTNAME=THROUGHPUT
----------------------

binary table with columns:

- ``wavelength[nwave]``   Wavelength in Angstroms, or `loglam[nwave]`  for log10(Angstroms)
- ``extinction[nwave]``   Atmospheric extinction in mags/airmass;
  Atm throughput = 10^(0.4 * extinction * airmass)
- ``throughput[nwave]``   Telescope, fibers, spectrographs, and CCDs
- ``fiberinput[nwave]``   Geometric loss at input of fiber for point source.

Required keywords in same HDU ``EXTNAME=THROUGHPUT`` (not HDU 0)

- ``EXPTIME``:  Standard exposure time in seconds
- ``EFFAREA``:  Effective area of mirror, including obscuration effects, in cm^2
- ``FIBERDIA``: Fiber diameter in arcsec

*Note*:
The throughtput format currently does not support per-fiber throughput.

*Note*:
``fiberinput`` in general depends upon source size, seeing, and fiber position
misalignments.  The default is for for perfectly centered point sources.
The format currently does not support extented objects (the user can mimic
this though the input spectra provided to Specter.)

The various forms of throughput affect different types of sources::

                 OBJECT   SKY    CALIB
    EXTINCTION    yes     yes     no
    FIBERINPUT    yes     no      no
    THROUGHPUT    yes     yes     yes

CALIB = calibration lamps or laser systems mounted with the telescope,
*i.e.* not experiencing sky absorption.
