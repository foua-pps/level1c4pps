#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c21529.ad.smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Unit tests for angle computations."""

from level1c4pps import make_azidiff_angle, update_angle_attributes

import datetime as dt
import numpy as np
import xarray as xr
try:
    from unittest import mock
except ImportError:
    import mock
import sys
if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

SAT_AZ = np.ma.array([[48.0,  56.0, 64.0,  72.0],
                      [80.0,  88.0, 96.0, 104.0],
                      [-80.0,  -88.0, -96.0, -104.0],
                      [-180.0,  -188.0, -196.0, -204.0],
                      [261.88999414630234,  266.1, 271.1, 276.1]], mask=False)
XSAT_AZ = xr.DataArray([[48.0,  56.0, 64.0,  72.0],
                        [80.0,  88.0, 96.0, 104.0],
                        [-80.0,  -88.0, -96.0, -104.0],
                        [-180.0,  -188.0, -196.0, -204.0],
                        [261.88999414630234, 266.1, 271.1, 276.1]])
SUN_AZ = np.ma.array([[148.0,  156.0, 164.0,  172.0],
                      [180.0,  188.0, 196.0, 204.0],
                      [180.0,  188.0, 196.0, 204.0],
                      [185.0,  193.0, 201.0, 209.0],
                      [196.77999560162425, 201.1, 206.1, 211.1]], mask=False)
XSUN_AZ = xr.DataArray([[148.0,  156.0, 164.0,  172.0],
                        [180.0,  188.0, 196.0, 204.0],
                        [180.0,  188.0, 196.0, 204.0],
                        [185.0,  193.0, 201.0, 209.0],
                        [196.77999560162425, 201.1, 206.1, 211.1]])

RES = np.ma.array([[100., 100., 100., 100.],
                   [100., 100., 100., 100.],
                   [100.,  84.,  68.,  52.],
                   [5.,  21.,  37.,  53.],
                   [65.10999854, 65., 65., 65.]],
                  mask=False)

XRES = xr.DataArray([[100., 100., 100., 100.],
                     [100., 100., 100., 100.],
                     [100.,  84.,  68.,  52.],
                     [5.,  21.,  37.,  53.],
                     [65.10999854, 65., 65., 65.]])


class TestAzimuthDifferenceAngles(unittest.TestCase):
    """Class for testing the derivation of azimuth difference angles."""

    def test_make_azidiff_angle(self):
        """Test calculating the azimuth difference angles."""
        daz = make_azidiff_angle(SAT_AZ, SUN_AZ)
        np.testing.assert_almost_equal(daz, RES)

        # Xarray data:
        xdaz = make_azidiff_angle(XSAT_AZ, XSUN_AZ)
        xr.testing.assert_allclose(xdaz, XRES, rtol=0.00001)


class TestUpdateAnglesAttribute(unittest.TestCase):
    """Test setting of attributes for angles."""

    def test_update_angle_attributes(self):
        """Test setting of attributes for angles."""
        band = mock.MagicMock(attrs={'start_time':  dt.datetime(2009, 7, 1, 12, 15),
                                     'end_time':  dt.datetime(2009, 7, 1, 12, 16)})

        class AngleObj(object):
            def __init__(self):
                self.attrs = {'area': 'xx'}
                self.coords = {}
        angle_obj = AngleObj()
        setattr(angle_obj, 'attrs', {'area': 'xx'})
        setattr(angle_obj, 'coords', {})
        angle_dict = {'satzenith': AngleObj(),
                      'sunzenith': AngleObj(),
                      'azimuthdiff': AngleObj()}
        update_angle_attributes(angle_dict, band)
        self.assertEqual(angle_dict['satzenith'].attrs['name'], 'satzenith')
        self.assertEqual(angle_dict['satzenith'].attrs['id_tag'], 'satzenith')
        np.testing.assert_array_equal(angle_dict['satzenith'].attrs['valid_range'],
                                      np.array([0, 9000], dtype='int16'))
        self.assertNotIn('area', angle_dict['satzenith'].attrs.keys())
        self.assertIn('time', angle_dict['satzenith'].coords.keys())
        self.assertIn('long_name', angle_dict['satzenith'].attrs.keys())
        self.assertIn('standard_name', angle_dict['satzenith'].attrs.keys())


def suite():
    """Create the test suite for test_atm_correction_ir."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestAzimuthDifferenceAngles))
    mysuite.addTest(loader.loadTestsFromTestCase(TestUpdateAnglesAttribute))
    return mysuite
