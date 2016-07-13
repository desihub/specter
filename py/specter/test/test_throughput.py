#!/usr/bin/env python

"""
Test specter throughput file format
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import numpy as np
import unittest
from pkg_resources import resource_filename
from ..throughput import load_throughput

class TestThroughput(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.thru = load_throughput(resource_filename('specter.test', 't/throughput.fits'))
        cls.w = np.arange(5000, 9000, 1)

    def setUp(self):
        self.flux = np.random.uniform(1,2, size=self.w.shape) * 1e-17
        self.photflux = np.random.uniform(0,1, size=self.w.shape)

    def test_area(self):
        self.assertTrue(self.thru.area > 0)

    def test_exptime(self):
        self.assertTrue(self.thru.exptime > 0)

    def test_fiberdia(self):
        self.assertTrue(self.thru.fiberdia > 0)

    def test_atmthru(self):
        t = self.thru.atmospheric_throughput(self.w)
        self.assertTrue(np.all( (0.0<=t) & (t<= 1.0)))

    def test_atmext(self):
        ext = self.thru.extinction(self.w)
        self.assertTrue(np.all( ext > 0 ))

    def test_fiberin(self):
        t = self.thru.fiberinput_throughput(self.w)
        self.assertTrue(np.all( (0.0<=t) & (t<= 1.0)))

        #- Use a different wavelength grid
        t = self.thru.fiberinput_throughput(self.w[0::2])
        self.assertTrue(np.all( (0.0<=t) & (t<= 1.0)))
        self.assertTrue(len(t) == len(self.w)//2)

        #- Should even work with no wavelength grid
        t = self.thru.fiberinput_throughput()
        self.assertTrue(np.all( (0.0<=t) & (t<= 1.0)))

    def test_fiberin_objtype(self):
        tstar = self.thru.fiberinput_throughput(self.w, objtype='STAR')
        telg = self.thru.fiberinput_throughput(self.w, objtype='ELG')
        tsky = self.thru.fiberinput_throughput(self.w, objtype='SKY')

        self.assertTrue(np.all( (0.0<=tstar) & (tstar<= 1.0)))
        self.assertTrue(np.all( (0.0<=telg) & (telg<= 1.0)))
        self.assertTrue(np.allclose(tsky, 1.0))

        self.assertTrue(np.all(tsky >= tstar))
        self.assertTrue(np.all(tstar >= telg))

    def test_hardware(self):
        t = self.thru.hardware_throughput(self.w)
        self.assertTrue(np.all( (0.0<=t) & (t<= 1.0)))

    def test_calibthru(self):
        t1 = self.thru(self.w, objtype='CALIB', airmass=1.0)
        t2 = self.thru(self.w, objtype='CALIB', airmass=2.0)
        self.assertTrue( np.all(t1 == t2) )
        self.assertTrue( len(t1) == len(t2) )
        self.assertTrue( len(t1) == len(self.w) )
        self.assertTrue( np.any(t2>0.0) )

    def test_skythru(self):
        t1 = self.thru(self.w, objtype='SKY', airmass=1.0)
        t2 = self.thru(self.w, objtype='SKY', airmass=2.0)
        self.assertTrue( np.all( (t1 > t2) | (t2 == 0.0)) )
        self.assertTrue( len(t1) == len(t2) )
        self.assertTrue( len(t1) == len(self.w) )
        self.assertTrue( np.any(t2>0.0) )

    def test_starthru(self):
        t1 = self.thru(self.w, objtype='STAR', airmass=1.0)
        t2 = self.thru(self.w, objtype='STAR', airmass=2.0)
        self.assertTrue( np.all( (t1 > t2) | (t2 == 0.0)) )
        self.assertTrue( len(t1) == len(t2) )
        self.assertTrue( len(t1) == len(self.w) )
        self.assertTrue( np.any(t2>0.0) )

    def test_thru(self):
        t1 = self.thru(self.w, objtype='CALIB', airmass=1.1)
        t2 = self.thru(self.w, objtype='SKY', airmass=1.1)
        t3 = self.thru(self.w, objtype='STAR', airmass=1.1)

        self.assertTrue( np.all( (t1>t2) | (t1==0.0) ) )
        self.assertTrue( np.all( (t2>t3) | (t2==0.0) ) )
        self.assertTrue( np.any(t3>0.0) )

    def test_fluxunits(self):
        units = [
            "erg/s/cm^2",
            "erg/s/cm^2/A",
            "erg/s/cm^2/A/arcsec^2",
            "erg/s/cm^2/arcsec^2",
        ]
        for u in units:
            p = self.thru.photons(self.w, self.flux, units=u, objtype='STAR')
            self.assertTrue( np.any(p>0) and np.all(p>=0) )

    def test_scaledunits(self):
        scale = 1e-16
        units = [
            "erg/s/cm^2",
            "erg/s/cm^2/A",
            "erg/s/cm^2/A/arcsec^2",
            "erg/s/cm^2/arcsec^2",
        ]
        for u in units:
            scaled_units = str(scale) + " " + u
            p0 = self.thru.photons(self.w, self.flux, units=u, objtype='STAR')
            p1 = self.thru.photons(self.w, self.flux/scale, units=scaled_units, objtype='STAR')
            ii = (p0 != 0.0)
            dp = np.abs( (p0-p1)[ii]/p0[ii] )
            self.assertTrue( np.max(dp) < 1e-14 )  #- allow for roundoff error

    def test_photunits(self):
        units = [
            "photon",
            "photon/A",
            ### "photon/A/arcsec^2",
        ]
        for u in units:
            p = self.thru.photons(self.w, self.photflux, units=u, objtype='STAR')
            self.assertTrue( np.any(p>0) and np.all(p>=0) )

    def test_calibphot(self):
        p1 = self.thru.photons(self.w, self.flux, objtype='CALIB', airmass=1.0)
        p2 = self.thru.photons(self.w, self.flux, objtype='CALIB', airmass=2.0)
        self.assertTrue( np.any(p1>0) and np.all(p1==p2) )

    def test_skyphot(self):
        p1 = self.thru.photons(self.w, self.flux, objtype='SKY', airmass=1.0)
        p2 = self.thru.photons(self.w, self.flux, objtype='SKY', airmass=2.0)
        self.assertTrue( np.any(p1>0) )
        self.assertTrue( np.all( (p1>p2) | (p2==0) ) )

    def test_objphot(self):
        p1 = self.thru.photons(self.w, self.flux, objtype='STAR', airmass=1.0)
        p2 = self.thru.photons(self.w, self.flux, objtype='STAR', airmass=2.0)
        self.assertTrue( np.any(p1>0) )
        self.assertTrue( np.all( (p1>p2) | (p2==0) ) )

    def test_multiphot(self):
        p1 = self.thru.photons(self.w, self.flux, objtype='CALIB', airmass=1.0)
        p2 = self.thru.photons(self.w, self.flux, objtype='SKY', airmass=1.0)
        p3 = self.thru.photons(self.w, self.flux, objtype='STAR', airmass=1.0)

        self.assertTrue( np.all( (p1>p2) | (p1==0.0) ) )
        self.assertTrue( np.all( (p2>p3) | (p2==0.0) ) )
        self.assertTrue( np.any(p3>0.0) )

    def test_apply_throughput(self):
        f1 = self.thru.apply_throughput(self.w, self.flux, objtype='CALIB')
        f2 = self.thru.apply_throughput(self.w, self.flux, objtype='SKY')
        f3 = self.thru.apply_throughput(self.w, self.flux, objtype='STAR')

        self.assertTrue( np.all( (f1>f2) | (f1==0.0) ) )
        self.assertTrue( np.all( (f2>f3) | (f2==0.0) ) )
        self.assertTrue( np.any(f3>0.0) )

        self.assertTrue( np.all(f1 <= self.flux) )
        self.assertTrue( np.all(f2 <= self.flux) )
        self.assertTrue( np.all(f3 <= self.flux) )

        self.assertTrue( np.all(f2 <= f1) )
        self.assertTrue( np.all(f3 <= f2) )

    def test_apply_throughput_multi(self):
        flux = np.array([self.flux, self.flux, self.flux])
        objtype = ['CALIB', 'SKY', 'STAR']
        fx = self.thru.apply_throughput(self.w, flux, objtype=objtype)
        self.assertEqual(fx.shape, flux.shape)
        self.assertTrue(np.all(fx <= self.flux))
        self.assertTrue( np.all(fx[1] <= fx[0]) )
        self.assertTrue( np.all(fx[2] <= fx[1]) )

    def test_wavemin_wavemax(self):
        wavemin = self.thru.wavemin
        wavemax = self.thru.wavemax
        self.assertLess(wavemin, wavemax)

    #- This test currently works, but some pixelated rebinning models can
    #- involve negative ringing which would cause this test to fail.
    def test_low_flux(self):
        flux = np.random.uniform(0, 1, len(self.w))
        f1 = self.thru.apply_throughput(self.w, flux, objtype='CALIB')
        f2 = self.thru.apply_throughput(self.w, flux, objtype='SKY')
        f3 = self.thru.apply_throughput(self.w, flux, objtype='STAR')

        self.assertTrue( np.all(f1 <= flux) )
        self.assertTrue( np.all(f2 <= flux) )
        self.assertTrue( np.all(f3 <= flux) )

        self.assertTrue( np.all(f2 <= f1) )
        self.assertTrue( np.all(f3 <= f2) )

if __name__ == '__main__':
    unittest.main()
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestThroughputIO)
    # unittest.TextTestRunner(verbosity=2).run(suite)
