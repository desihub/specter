from __future__ import absolute_import, division, print_function, unicode_literals

import unittest
import numpy as np
from specter.util import pixspline

class TestPixSpline(unittest.TestCase):
    """
    Test functions within specter.util
    """

    def test_cen2bounds(self):
        #- Integers
        x = np.arange(10)
        xedge = pixspline.cen2bound(x)
        dxedge = np.diff(xedge)
        self.assertTrue(np.all(dxedge > 0))
        self.assertTrue(np.all(dxedge == dxedge[0]))

        #- Floats
        x = np.arange(10.0)
        xedge = pixspline.cen2bound(x)
        dxedge = np.diff(xedge)
        self.assertTrue(len(xedge) == len(x)+1)
        self.assertTrue(np.all(dxedge > 0))
        self.assertTrue(np.all(dxedge == dxedge[0]))

    def test_pixspline(self):
        nx = 100
        x = np.arange(nx)
        y = np.random.uniform(0,1, size=nx)
        edges = pixspline.cen2bound(x)

        #- PixelSpline with edges
        sp = pixspline.PixelSpline(edges, y)
        dy = sp.resample(edges) - y
        self.assertTrue(np.max(np.abs(dy)) < 1e-12)

        #- PixelSpline with centers
        sp = pixspline.PixelSpline(x, y)
        dy = sp.resample(edges) - y
        self.assertTrue(np.max(np.abs(dy)) < 1e-12)

        #- Test with vectors vs scalars
        yy = sp.point_evaluate(x)
        for i in range(len(yy)):
            self.assertTrue(sp.point_evaluate(x[i]) == yy[i])

        #- Rebinning
        y2 = y.reshape((nx//2,2)).mean(1)
        dy2 = sp.resample(edges[0::2]) - y2
        self.assertTrue(np.max(np.abs(dy2)) < 1e-12)

        y5 = y.reshape((nx//5,5)).mean(1)
        dy5 = sp.resample(edges[0::5]) - y5
        self.assertTrue(np.max(np.abs(dy5)) < 1e-12)

        #- Sanity check that other functions don't crash
        blat = sp.find_extrema()
        blat = sp.point_evaluate(x)
        blat = sp(x)

if __name__ == '__main__':
    unittest.main()
