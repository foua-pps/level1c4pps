#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 level1c4pps developers
#
# This file is part of level1c4pps
#
# level1c4pps is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# level1c4pps is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with level1c4pps.  If not, see <http://www.gnu.org/licenses/>.

"""Unit tests for the fci2pps_lib module."""
import datetime as dt
import sys
import unittest

import numpy as np
import pytest
import xarray as xr
from pyresample import SwathDefinition
from satpy.dataset.dataid import WavelengthRange

try:
    from unittest import mock
except ImportError:
    import mock
sys.modules["hdf5plugin"] = mock.MagicMock()
import level1c4pps.fci2pps_lib as fci2pps  # modify hdf5plugin first # noqa: E402


def get_fake_scene(start_time=dt.datetime(2000, 1, 1, 0)):
    """Create fake scene."""
    from satpy import Scene
    scene = Scene()
    end_time = dt.datetime(2000, 1, 1, 0, 1)
    lons = xr.DataArray(
        [[75.0, 76.0],
         [77.0, 78.0]],
        dims=('y', 'x'))
    lats = xr.DataArray(
        [[35.0, 36.0],
         [37.0, 38.0]],
        dims=('y', 'x'))
    scene['vis_06'] = xr.DataArray(
        [[1.0, 2.0],
         [3.0, 4.0]],
        dims=('y', 'x'),
        attrs={'calibration': 'reflectance',
               'sun_earth_distance_correction_applied': True,
               'start_time': start_time,
               'end_time': end_time,
               'wavelength': WavelengthRange(0.56, 0.635, 0.71)}
    )
    scene['ir_105'] = xr.DataArray(
        [[5, 6],
         [7, 8]],
        dims=('y', 'x'),
        attrs={'calibration': 'brightness_temperature',
               'start_time': start_time,
               'end_time': end_time,
               'wavelength': WavelengthRange(9.8, 10.8, 11.8),
               'area': SwathDefinition(lons, lats)}
    )
    scene['ir_105_time'] = xr.DataArray(
        [[5.0, 6.0],
         [7.2, 8.0]],
        dims=('y', 'x'),
        attrs={}

    )
    scene.attrs['sensor'] = {'fci'}
    scene.attrs['platform_name'] = 'Meteosat-12'
    return scene


class TestFCI2PPS(unittest.TestCase):
    """Test for FCI converter."""

    def test_fix_time(self):
        """Test loading and calibrating the data."""
        scene = get_fake_scene()
        fci2pps.fix_time(scene)
        self.assertEqual(scene['ir_105_time'][0, 0], np.datetime64('2000-01-01T00:00:05'))
        scene = get_fake_scene(start_time=dt.datetime(2000, 1, 1, 0, 0, 59))
        with self.assertRaises(ValueError):
            fci2pps.fix_time(scene)

    @mock.patch('satpy.utils.get_satpos')
    def test_add_angles_and_latlon(self, my_get_satpos):
        """Test adding angles and lat/lon."""
        my_get_satpos.return_value = (-0.3, 0.0, 35786 * 1000)
        scene = get_fake_scene()
        fci2pps.add_angles_and_latlon(scene)
        np.testing.assert_allclose(scene["satzenith"].values[0, 0], 86.6785348)

    def test_resample(self):
        """Test resampling."""
        scene = get_fake_scene()
        out1 = fci2pps.resample_data(scene, ["ir_105"], "coarse")
        out2 = fci2pps.resample_data(scene, ["ir_105"], "fine")
        out3 = fci2pps.resample_data(scene, ["ir_105"], "msg_seviri_fes_3km")
        self.assertEqual(out3["ir_105"].shape[0], 3712)
        self.assertEqual(out1["ir_105"].shape[0], 2)
        self.assertEqual(out2["ir_105"].shape[0], 2)

    def test_attrs(self):
        """Test setting of attributes."""
        scene = get_fake_scene()
        fci2pps.set_header_and_band_attrs(scene)
        self.assertEqual(scene["vis_06"].attrs["description"], "FCI VIS_06")
        self.assertEqual(scene["vis_06"].attrs["name"], "image1")
        self.assertEqual(scene["vis_06"].attrs["sun_zenith_angle_correction_applied"], "True")
        self.assertEqual(scene["vis_06"].attrs["valid_range"][1], 20000)


def suite():
    """Create the test suite for test_fci2pps."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestFCI2PPS))
    return mysuite
