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
import os

import numpy as np
import xarray as xr
from pyresample import SwathDefinition
from satpy.dataset.dataid import WavelengthRange

from unittest import mock
sys.modules["hdf5plugin"] = mock.MagicMock()
import level1c4pps.fci2pps_lib as fci2pps  # modify hdf5plugin first # noqa: E402


def get_fake_scene(start_time=dt.datetime(2000, 1, 1, 0)):
    """Create fake scene."""
    from satpy import Scene
    scene = Scene()
    end_time = dt.datetime(2000, 1, 1, 0, 1)
    lons = xr.DataArray(
        74 + np.ones((6000, 6000)),
        dims=('y', 'x'))
    lats = xr.DataArray(
        34 + np.ones((6000, 6000)),
        dims=('y', 'x'))
    scene['vis_06'] = xr.DataArray(
        np.ones((6000, 6000)),
        dims=('y', 'x'),
        attrs={'calibration': 'reflectance',
               'sun_earth_distance_correction_applied': True,
               'start_time': start_time,
               'end_time': end_time,
               'wavelength': WavelengthRange(0.56, 0.635, 0.71),
               'area': SwathDefinition(lons, lats)}
    )
    scene['vis_09'] = xr.DataArray(
        np.ones((3000, 3000)),
        dims=('y', 'x'),
        attrs={'calibration': 'reflectance',
               'sun_earth_distance_correction_applied': True,
               'start_time': start_time,
               'end_time': end_time,
               'wavelength': WavelengthRange(0.56, 0.635, 0.71),
               'area': SwathDefinition(lons[::2, ::2], lats[::2, ::2])}
    )
    scene['ir_105_time'] = xr.DataArray(
        np.ones((6000, 6000)),
        dims=('y', 'x'),
        attrs={"name": 'ir_105_time'}
    )
    scene["longitude"] = lons
    scene["latitude"] = lats
    scene.attrs["end_time"] = end_time
    scene.attrs["start_time"] = start_time
    scene['ir_105'] = xr.DataArray(
        4.0 * np.ones((6000,6000)),
        dims=('y', 'x'),
        attrs={'calibration': 'brightness_temperature',
               'start_time': start_time,
               'end_time': end_time,
               'wavelength': WavelengthRange(9.8, 10.8, 11.8),
               'area': SwathDefinition(lons, lats)}
    )
    scene['ir_105_time'] = xr.DataArray(
        4.0 + np.ones((6000,6000)),
        dims=('y', 'x'),
        attrs={}
    )
    scene.attrs['sensor'] = {'fci'}
    scene.attrs['platform_name'] = 'Meteosat-12'
    scene.load = mock.MagicMock
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
        self.assertEqual(out1["ir_105"].shape[0], 3000)
        self.assertEqual(out2["ir_105"].shape[0], 6000)
        with self.assertRaises(ValueError):
            fci2pps.resample_data(scene, ["ir_105"], "bad_arg")

    def test_attrs(self):
        """Test setting of attributes."""
        scene = get_fake_scene()
        fci2pps.set_header_and_band_attrs(scene)
        self.assertEqual(scene["vis_06"].attrs["description"], "FCI VIS_06")
        self.assertEqual(scene["vis_06"].attrs["name"], "image1")
        self.assertEqual(scene["vis_06"].attrs["sun_zenith_angle_correction_applied"], "True")
        self.assertEqual(scene["vis_06"].attrs["valid_range"][1], 20000)

    @mock.patch("level1c4pps.fci2pps_lib.Scene")
    @mock.patch("level1c4pps.fci2pps_lib.add_angles_and_latlon")
    def test_process_one_scene(self, mock_lonlats, mock_scene):
        """Test to set process_one_scene."""
        import level1c4pps.fci2pps_lib as fci2pps
        scene = get_fake_scene()
        del scene["vis_09"]
        mock_scene.return_value = scene
        scene.resample_data = mock.MagicMock(return_value=scene)
        filename = fci2pps.process_one_scene("dummy", out_path='./level1c4pps/tests/')
        self.assertEqual(os.path.basename(filename),
                         "S_NWC_fci_meteosat12_00000_20000101T0000000Z_20000101T0001000Z.nc")
