# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
specter.test.specter_test_suite
===============================

Used to initialize the unit test framework via ``python setup.py test``.
"""
#
from __future__ import absolute_import, division, print_function, unicode_literals
#
import unittest
#
#- This is factored out separately from runtests() so that it can be used by
#- python setup.py test
def specter_test_suite():
    """Returns unittest.TestSuite of specter tests"""
    from os.path import dirname
    specter_dir = dirname(dirname(__file__))
    # print(specter_dir)
    return unittest.defaultTestLoader.discover(specter_dir,
        top_level_dir=dirname(specter_dir))

def runtests():
    """Run all tests in specter.test.test_*.
    """
    #- Load all TestCase classes from specter/test/test_*.py
    tests = specter_test_suite()
    #- Run them
    unittest.TextTestRunner(verbosity=2).run(tests)

if __name__ == "__main__":
    runtests()
