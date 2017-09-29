=====================
specter Release Notes
=====================

0.8.1 (unreleased)
------------------

* No changes yet.

0.8.0 (2017-09-29)
------------------

* Added subbundle divide-and-conquer extractions for ~2x speedup

0.7.0 (2017-03-02)
------------------

* Update template Module file to reflect DESI+Anaconda infrastructure.
* Enable projecting photons onto multiple images simultaneously
* Fix GaussHermite PSF spot size and centering bugs
* New PSF function ._value to evaluate non-pixel-integrated PSF values

0.6.0 (2016-08-16)
------------------

PR #40:
* Added full_output option to ex2d to get model image and metrics based upon
  goodness of fit
* PSFs can specify their model error with PSFERR header keyword; default 0.01

0.5.0 (2016-05-23)
------------------

* Move data files into Python package so pip can install the data files.
* Load test files in class methods to hopefully speed up tests.
* Improve Travis test support to latest standards.
* Added a documentation page for the specter API.

0.4.1 (2016-03-10)
------------------

* Bug fixes for small PSFs, and fixes of the fixes
* This is a release candidate for DESI Spectro Pipeline 2016a

0.4 (2016-03-03)
----------------

* refactored bin/exspec to move most functionality into specter.extract.ex2d
  API change to ex2d() to use specmin,nspec instead of
  specrange=(specmin,specmax)
* removed desiutil dependency

0.3 (2015-12-15)
----------------

* pip install support, among many changes.
* This version includes the desiutil infrastructure.  This will probably be
  removed in the future, but for now this is needed for installation support.

0.2.5 (2015-04-14)
------------------

* Includes cachedict bug fix and traceset.fit_traces utility function.

0.2.4 (2015-02-13)
------------------

* "robot overlords"
* use scipy.linalg instead of numpy.linalg

0.2.3 (2015-02-05)
------------------

* more linalg stability attempts
* ivar renaming typo

0.2.2 (2015-02-03)
------------------

* trim by percent of median not percentile

0.2.1 (2015-02-02)
------------------

* Added better (?) linear algebra conditioning; dump a debug file if the linear algebra fails.

0.2 (2015-02-02)
----------------

* GaussHermite vs. GaussHermite2 from dev branch

0.1.3 (2015-01-24)
------------------

* More robust when pixels are masked
* Adds a linear algebra robustness check for when pixels are masked or when asking for wavelengths that are entirely off the CCD.

0.1.2 (2015-01-07)
------------------

* Fixes a bug when asking for xyrange for wavelengths that are way off the CCD and the extrapolation has gone very bad.

0.1.1 (2015-01-06)
------------------

* Bug fix to xyrange when wavelengths are within a half a pixel of the CCD boundary.

0.1 (2014-12-29)
----------------

* Initial tag.
