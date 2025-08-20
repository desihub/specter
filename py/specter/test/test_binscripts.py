#!/usr/bin/env python

"""
Unit tests for executable scripts in specter/bin
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import sys
import unittest
import subprocess as sp
from shutil import rmtree
from tempfile import mkdtemp
from uuid import uuid4
from importlib.resources import files
import numpy as np
from astropy.io import fits
from specter.io import read_image, write_spectra


class TestBinScripts(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """This method runs before any other methods in this class.
        """
        cls.tmp_dir = mkdtemp()
        _base = uuid4().hex
        cls.imgfile1 = os.path.join(cls.tmp_dir, f'testimg1-{_base}.fits')
        cls.imgfile2 = os.path.join(cls.tmp_dir, f'testimg2-{_base}.fits')
        cls.imgfile_exspec = os.path.join(cls.tmp_dir, f'testimg-exspec-{_base}.fits')
        cls.specfile1 = os.path.join(cls.tmp_dir, f'testspec1-{_base}.fits')
        cls.specfile2 = os.path.join(cls.tmp_dir, f'testspec2-{_base}.fits')

        #- when running "python setup.py test", this file is run from different
        #- locations for python 2.7 vs. 3.5
        #- python 2.7: py/specter/test/test_binscripts.py
        #- python 3.5: build/lib/specter/test/test_binscripts.py

        #- python 2.7 location:
        cls.specter_dir = os.path.dirname( # top-level
            os.path.dirname( # py/
                os.path.dirname( # specter/
                    os.path.dirname(os.path.abspath(__file__)) # test/
                    )
                )
            )
        if not os.path.isdir(cls.specter_dir + '/bin'):
            #- python 3.x setup.py test location:
            cls.specter_dir = os.path.dirname( # top-level
                os.path.dirname( # build/
                    os.path.dirname( # lib/
                        os.path.dirname( # specter/
                            os.path.dirname(os.path.abspath(__file__)) # test/
                            )
                        )
                    )
                )
        #- last attempt
        if not os.path.isdir(cls.specter_dir + '/bin'):
            cls.specter_dir = os.path.join(os.getcwd(), 'bin')

        if not os.path.isdir(cls.specter_dir + '/bin'):
            raise RuntimeError('Unable to auto-locate specter/bin from {}'.format(__file__))

        cls.executable = sys.executable
        cls.sky_file = str(files('specter').joinpath('data', 'sky-uves.fits'))
        cls.monospot_file = str(files('specter').joinpath('test', 't', 'psf-monospot.fits'))
        cls.throughput_file = str(files('specter').joinpath('test', 't', 'throughput.fits'))

        #- Add this package to PYTHONPATH so that binscripts can find it
        try:
            cls.origPath = os.environ['PYTHONPATH']
            os.environ['PYTHONPATH'] = os.path.join(cls.specter_dir,'py') + ':' + cls.origPath
        except KeyError:
            cls.origPath = None
            os.environ['PYTHONPATH'] = os.path.join(cls.specter_dir,'py')

        # Generate input file needed for exspec tests.
        cmd = cls.specter_command(cls.imgfile_exspec, extra=True, noise=True)
        proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = proc.communicate()
        cls.extra_returned = (proc.returncode, out.decode('utf-8'), err.decode('utf-8'))

    @classmethod
    def tearDownClass(cls):
        """This method runs after all other test methods in this class.
        """
        if cls.origPath is None:
            del os.environ['PYTHONPATH']
        else:
            os.environ['PYTHONPATH'] = cls.origPath

        rmtree(cls.tmp_dir)

    def setUp(self):
        """This method runs before every individual test method in this class.
        """
        pass

    def tearDown(self):
        """This method runs after every individual test method in this class.
        """
        pass

    @classmethod
    def specter_command(cls, imagefile, extra=False, noise=False, trimxy=False):
        """Generate a specter command.

        Parameters
        ----------
        imagefile : :class:`str`
            2d image file name used as output.
        extra : :class:`bool`, optional
            If ``True``, add ``--extra`` to the command.
        noise : :class:`bool`, optional
            If ``True``, add ``--noise`` to the command.
        trimxy : :class:`bool`, optional
            If ``True``, add ``--trimxy`` to the command.

        Returns
        -------
        :class:`list`
            A list suitable for passing to :class:`subprocess.Popen`.
        """
        cmd = [cls.executable,
               os.path.join(cls.specter_dir, 'bin', 'specter'),
               '-i', cls.sky_file,  # Input spectra, in this case sky lines.
               '-p', cls.monospot_file,  # Input PSF
               '-t', cls.throughput_file,  # Throughput file
               '-o', imagefile,  # Output file
               '-w', '7500,7620',  # Wavelength range
               '--specmin', '0',  # Starting index for simulated spectra
               '--nspec', '2',  # Number of spectra to simulate
               '--exptime', '1500']  # Exposure time
        if extra:
            cmd += ['--extra']  # Write extra HDUs to output.
        if noise:
            cmd += ['--noise']  # Add noise.
        if trimxy:
            cmd += ['--trimxy']  # Trim output image.
        return cmd

    def exspec_command(self, imagefile, specfile, dwave, specmin=0, nspec=2):
        """Generate a exspec command.

        Parameters
        ----------
        imagefile : :class:`str`
            "Raw" image file used as input.
        specfile : :class:`str`
            "Extracted" spectrum file used as output.
        dwave : :class:`float`
            Interval between wavelength bins.
        specmin : :class:`int`, optional
            First spectrum to extract.
        nwave : :class:`int`, optional
            Number of spectra to extract.

        Returns
        -------
        :class:`list`
            A list suitable for passing to :class:`subprocess.Popen`.
        """
        cmd = [self.executable,
               os.path.join(self.specter_dir, 'bin', 'exspec'),
               '-i', imagefile,  # Input image
               '-p', self.monospot_file,  # Input PSF
               '-o', specfile,  # output extracted spectra
               '-w', f'7500,7620,{dwave}',  # wavemin,wavemax,dw
               '--specmin', str(specmin),  # first spectrum to extract
               '--nspec', str(nspec)]  # Number of spectra to extract
        return cmd

    def test_specter_command(self):
        """Test the outputs of a specter command.
        """
        cmd = self.specter_command(self.imgfile1, noise=True)
        # print(cmd)
        proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = proc.communicate()
        self.assertEqual(proc.returncode, 0, 'Error code {} != 0'.format(proc.returncode))
        self.assertTrue(os.path.exists(self.imgfile1))

        with fits.open(self.imgfile1) as fx:
            self.assertIn('CCDIMAGE', fx)
            self.assertIn('IVAR', fx)

        #- Test the I/O routines while we have the file handy
        image, ivar, hdr = read_image(self.imgfile1)
        self.assertEqual(image.shape, ivar.shape)
        self.assertTrue(image.dtype.isnative)
        self.assertTrue(ivar.dtype.isnative)

        os.remove(self.imgfile1)

    def test_specter_command_extra(self):
        """Test specter ... --extra.
        """
        self.assertEqual(self.extra_returned[0], 0, 'Error code {} != 0'.format(self.extra_returned[0]))
        self.assertTrue(os.path.exists(self.imgfile_exspec))
        with fits.open(self.imgfile_exspec) as fx:
            self.assertIn('PHOTONS', fx)
            self.assertIn('XYWAVE', fx)

    def test_specter_command_cores(self):
        """Compare specter results running on 1 or 2 cores.
        """
        cmd1 = self.specter_command(self.imgfile1, trimxy=True) + ['--numcores', '1']
        cmd2 = self.specter_command(self.imgfile2, trimxy=True) + ['--numcores', '2']

        if os.path.exists(self.imgfile1):
            os.remove(self.imgfile1)
        if os.path.exists(self.imgfile2):
            os.remove(self.imgfile2)

        proc = sp.Popen(cmd1, stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = proc.communicate()
        self.assertEqual(proc.returncode, 0, 'Error code {} != 0'.format(proc.returncode))
        self.assertTrue(os.path.exists(self.imgfile1))

        proc = sp.Popen(cmd2, stdout=sp.PIPE, stderr=sp.PIPE)
        out, err = proc.communicate()
        self.assertEqual(proc.returncode, 0, 'Error code {} != 0'.format(proc.returncode))
        self.assertTrue(os.path.exists(self.imgfile2))

        img1 = fits.getdata(self.imgfile1)
        img2 = fits.getdata(self.imgfile2)

        self.assertTrue(np.allclose(img1, img2))

    def test_exspec(self):
        """Test the exspec command.

        This test requires a file that should be generated by setUpClass().
        """
        check_specfile = True
        for dwave, specmin, nspec in [(1.0, 0, 2), (2.0, 0, 2), (1.0, 498, 2)]:
            cmd = self.exspec_command(self.imgfile_exspec, self.specfile1,
                                      dwave=dwave, specmin=specmin, nspec=nspec)
            proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
            out, err = proc.communicate()
            self.assertEqual(proc.returncode, 0,
                             'Error code {} != 0 with dwave={}, specmin={}, nspec={}.'.format(proc.returncode, dwave, specmin, nspec))
            self.assertTrue(os.path.exists(self.specfile1))

            if check_specfile:
                with fits.open(self.specfile1) as fx:
                    # print(fx.info())
                    self.assertIn('FLUX', fx)
                    self.assertIn('IVAR', fx)
                    self.assertIn('WAVELENGTH', fx)
                    self.assertIn('RESOLUTION', fx)

                    #- this is covered in the exspec binscript, but not yet visible to
                    #- coverage tools; try it here just for good measure
                    write_spectra(self.specfile2,
                        fx['WAVELENGTH'].data, fx['FLUX'].data,
                        fx['IVAR'].data, fx['RESOLUTION'].data, fx[0].header)
                check_specfile = False  # This check only needs to run once.
