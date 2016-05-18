from specter.test.test_psf import TestPixPSF, TestSpotPSF
from specter.test.test_specio import TestSpecIO
from specter.test.test_throughput import TestThroughput
from specter.test.test_util import TestUtil
from specter.test.test_extract import TestExtract
from specter.test.test_pixspline import TestPixSpline
from specter.test.test_binscripts import TestBinScripts

def test():
    """
    Run a suite of specter tests
    """
    import unittest

    load = unittest.defaultTestLoader.loadTestsFromTestCase
    tests = list()
    tests.append(load(TestPixPSF))
    tests.append(load(TestSpotPSF))
    tests.append(load(TestSpecIO))
    tests.append(load(TestThroughput))
    tests.append(load(TestUtil))
    tests.append(load(TestExtract))
    tests.append(load(TestPixSpline))
    tests.append(load(TestBinScripts))

    suite = unittest.TestSuite(tests)
    unittest.TextTestRunner(verbosity=2).run(suite)
