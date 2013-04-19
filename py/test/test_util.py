#!/usr/bin/env python

"""
Unit tests for PSF classes.
"""

import sys
import os
import numpy as N
from numpy.polynomial import legendre
import specter.util
import unittest

class TestUtil(unittest.TestCase):
    """
    Test functions within specter.util
    """
    
    def test_gaussX(self):
        #- Gaussian integration
        self.assertTrue(specter.util.gaussint(-100) == 0.0)
        self.assertTrue(specter.util.gaussint(0.0) == 0.5)
        self.assertTrue(specter.util.gaussint(+100) == 1.0)
        
        for x in (-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 10.0):
            self.assertTrue(specter.util.gaussint(x, mean=x, sigma=1.0) == 0.5)
            self.assertTrue(specter.util.gaussint(x, mean=x, sigma=2.0) == 0.5)
            if x>0:
                self.assertTrue(specter.util.gaussint(x, sigma=x) == 0.84134474606854293)
                self.assertTrue(specter.util.gaussint(-x, sigma=x) == 0.15865525393145707)
            

if __name__ == '__main__':
    unittest.main()           
