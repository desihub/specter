=====================
specter Release Notes
=====================

0.3.1 (unreleased)
------------------

* No changes yet.

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
