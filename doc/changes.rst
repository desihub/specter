==================
specter change log
==================

0.9.5 (unreleased)
------------------

* No changes yet

0.9.4 (2020-08-03)
------------------

* Update documentation configuration for ReadTheDocs (PR `#77`_).
* Fix last-bin bug in pixellated Gauss-Hermite integration (PR `#79`_).
* Astropy deprecation rename clobber -> overwrite (PR `#81`_).
* Add reproducibility unit tests (PR `#81`_).

.. _`#77`: https://github.com/desihub/specter/pull/77
.. _`#79`: https://github.com/desihub/specter/pull/79
.. _`#81`: https://github.com/desihub/specter/pull/81

0.9.3 (2020-04-16)
------------------

* Improve handling of heavily (or completely) masked inputs (PR #78).

0.9.2 (2020-04-07)
------------------

* Fix NaN flux coming from masked input pixels (PR #76).

0.9.1 (2018-11-07)
------------------

* Faster Pixelated Gauss-Hermite (pgh) for ~25-30% extraction speedup
  (PR #71 and #73).
* Faster xypix 10-50% extraction speedup (PR #74).
* Memory and operation order improvements for ~few percent speedup (PR #75).

0.9.0 (2018-09-26)
------------------

* Faster extractions by vectorizing and caching legval calls (PR #70).

0.8.7 (2018-07-26)
------------------

* Add custom `xsigma` and `ysigma` functions to GaussHermitePSF (PR #66).
* Don't use numba caching due to MPI race condition (PR #67).
* Small speed improvements (PR #68 and #69).

0.8.6 (2018-06-27)
------------------

* Added numba-ized legval for ~20% overall ex2d speedup (PR #61).
* Fixed tests (PR #62).
* Less regularization for ringing to lower bias (PR #63).

0.8.5 (2018-05-10)
------------------

* Allow user to override psferr in ex2d (PR #60)

0.8.4 (2018-03-29)
------------------

* np.outer replacement for 6% faster runtime (PR #58)

0.8.3 (2018-02-23)
------------------

* SpotGrid speedup (used by DESI pixsim); adds numba dependency (PR #56)

0.8.2 (2017-12-20)
------------------

* Don't require 2to3 during installations; fix license (PR #55)

0.8.1 (2017-10-25)
------------------

* Robust even if nsubbundles>bundlesize (PR #53)

0.8.0 (2017-09-29)
------------------

* Added subbundle divide-and-conquer extractions for ~2x speedup (PR #51)
* Added GaussHermite PSF format v3 (PR #52)

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
