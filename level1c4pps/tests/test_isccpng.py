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

"""Unit tests for the isccpng2pps_lib module."""

import datetime as dt
import unittest
import os
import numpy as np
from unittest import mock
from satpy import Scene
import xarray as xr

import level1c4pps.isccpng2pps_lib as isccpng2pps


class TestIsccpng2PPS(unittest.TestCase):
    """Test isccpng2pps_lib."""

    def setUp(self):
        """Create a test scene."""
        self.scene = Scene()
        scene_dict = {}
        grid_data = [[1.0, 2.0], [3.0, 4.0]]
        all_keys = ["refl_00_65um", 'temp_11_00um'] + isccpng2pps.GEOLOCATION_NAMES
        for key in all_keys:
            scene_dict[key] = xr.DataArray(grid_data,
                                           dims=('y', 'x'),
                                           attrs={'name': key,
                                                  'id_tag': key})
        scene_dict["refl_00_65um"].attrs = {'name': 'image0',
                                            'wavelength': [1, 2, 3, 'um'],
                                            'start_time': np.datetime64('2021-06-28T01:00:00.000000000+0100'),
                                            'end_time': np.datetime64('2021-06-28T01:01:00.000000000+0100'),
                                            'id_tag': 'ch_r06'}
        scene_dict['temp_11_00um'].attrs = {'name': 'image1',
                                            'id_tag': 'ch_tb11',
                                            'wavelength': [1, 2, 3, 'um'],
                                            'start_time': np.datetime64('2021-06-28T01:00:00.000000000+0100'),
                                            'end_time': np.datetime64('2021-06-28T01:01:00.000000000+0100'),
                                            'platform_name': '',
                                            'orbit_number': 99999}
        scene_dict["pixel_time"].coords["crs"] = ""
        for key in scene_dict:
            self.scene[key] = scene_dict[key]
        self.scene.attrs['sensor'] = ['isccpng']

    @mock.patch("level1c4pps.isccpng2pps_lib.load_data")
    def test_process_one_scene(self, mock_load):
        """Test to set process_one_scene."""
        import level1c4pps.isccpng2pps_lib as isccpng2pps
        mock_load.return_value = self.scene
        filename = isccpng2pps.process_one_scene("dummy", out_path='./level1c4pps/tests/', orbit_n='12345')
        self.assertEqual(os.path.basename(filename), "S_NWC_seviri__12345_20210628T0000000Z_20210628T0001000Z.nc")
