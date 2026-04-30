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

"""Unit tests for the metimage2pps_lib module."""

import datetime as dt
import unittest
import numpy as np
import xarray as xr
from satpy import Scene
from satpy.dataset.dataid import WavelengthRange
import level1c4pps.metimage2pps_lib as metimage2pps_lib


class TestMETimage2PPS(unittest.TestCase):
    """Test metimage2pps_lib."""

    def setUp(self):
        """Create a test scene."""
        metimage2pps_lib.BANDNAMES = ["vii_668",
                                      "vii_10690",
                                      "vii_12020"]
        scene = Scene()
        start_time = dt.datetime(2020, 1, 1, 12, 1)
        end_time = dt.datetime(2020, 1, 1, 12, 2)
        data = 200 + 100 * np.random.random((816, 3144))
        linear_48 = np.array(list(range(48)))
        linear_48 = linear_48 - np.mean(linear_48)
        scene["vii_668"] = xr.DataArray(
            data,
            dims=("y", "x"),
            attrs={"calibration": "reflectance",
                   "sun_earth_distance_correction_applied": False,
                   "start_time": start_time,
                   "wavelength": WavelengthRange(0.56, 0.635, 0.71)}
        )
        scene["vii_10690"] = xr.DataArray(
            data + np.tile(linear_48[:, np.newaxis], (17, 3144)),
            dims=("y", "x"),
            attrs={"calibration": "brightness_temperature",
                   "platform_name": "metopd",
                   "instrument": "metimage",
                   "start_time": start_time,
                   "end_time": end_time,
                   "wavelength": WavelengthRange(9.8, 10.8, 11.8)}
        )
        scene["vii_12020"] = xr.DataArray(
            data + 2 * np.tile(linear_48[:, np.newaxis], (17, 3144)),
            dims=("y", "x"),
            attrs={"calibration": "brightness_temperature",
                   "start_time": start_time,
                   "wavelength": WavelengthRange(9.8, 10.8, 11.8)}
        )
        self.scene = scene
        self.data = data

    def test_destripe(self):
        """Test destriping for METimage."""
        metimage2pps_lib.destripe(self.scene, "vii_12020", num=48)
        metimage2pps_lib.destripe(self.scene, "vii_10690", num=48)
        np.testing.assert_allclose(self.scene["vii_10690"], self.data, atol=0.5)
        np.testing.assert_allclose(self.scene["vii_12020"], self.data, atol=0.5)

    def test_set_header_and_band_attrs(self):
        """Test to set header_and_band_attrs."""
        metimage2pps_lib.set_header_and_band_attrs(self.scene)
        self.assertTrue(isinstance(self.scene.attrs["orbit_number"], int))
        self.assertEqual(self.scene["vii_668"].attrs["sun_zenith_angle_correction_applied"], "False")


def suite():
    """Create the test suite for test_metimage2pps."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestMETimage2PPS))

    return mysuite
