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

"""Unit tests for the gac2pps_lib module."""

import datetime as dt
import unittest
from unittest import mock
from satpy import Scene
import os
import xarray as xr
import level1c4pps.modis2pps_lib as modis2pps


class TestModis2PPS(unittest.TestCase):
    """Test modis2pps_lib."""

    def setUp(self):
        """Create a test scene."""
        self.scene = Scene()
        scene_dict = {}
        grid_data = [[1.0, 2.0], [3.0, 4.0]]
        all_keys = ['1', '31'] + modis2pps.GEOLOCATION_NAMES
        for key in all_keys:
            scene_dict[key] = xr.DataArray(grid_data,
                                           dims=('y', 'x'),
                                           attrs={'name': key,
                                                  'id_tag': key})
        scene_dict['1'].attrs = {'name': 'image0',
                                 'wavelength': [1, 2, 3, 'um'],
                                 'modifiers': ("sunz_correction", ),
                                 'id_tag': 'ch_r06'}
        scene_dict['31'].attrs = {'name': 'image1',
                                  'id_tag': 'ch_tb11',
                                  'platform': "Aqua",
                                  'number_of_scans': 15,
                                  'rows_per_scan': 16,
                                  'wavelength': [1, 2, 3, 'um'],
                                  'start_time': dt.datetime(2009, 7, 1, 12, 1, 0),
                                  'end_time': dt.datetime(2009, 7, 1, 12, 1, 0),
                                  'orbit_number': 99999}
        for key in scene_dict:
            self.scene[key] = scene_dict[key]
        self.scene.load = mock.MagicMock
        self.scene.attrs['sensor'] = ['modis']

    def test_compose_filename(self):
        """Test compose filename for MODIS."""
        start_time = dt.datetime(2009, 7, 1, 12, 15)
        end_time = dt.datetime(2009, 7, 1, 12, 30)
        scene = mock.MagicMock(attrs={'start_time': start_time,
                                      'end_time': end_time,
                                      'orbit_number': '99999',
                                      'platform': 'EOS-Aqua'})
        start_time = dt.datetime(2009, 7, 1, 12, 16)
        end_time = dt.datetime(2009, 7, 1, 12, 27)
        band = mock.MagicMock(attrs={'start_time': start_time,
                                     'end_time': end_time})
        fname_exp = '/out/path/S_NWC_modis_eos2_99999_20090701T1216000Z_20090701T1227000Z.nc'
        fname = modis2pps.compose_filename(scene, '/out/path', 'modis', band=band)
        self.assertEqual(fname, fname_exp)

    def test_set_header_and_band_attrs(self):
        """Test to set header_and_band_attrs."""
        modis2pps.set_header_and_band_attrs(self.scene, orbit_n='12345')
        self.assertTrue(isinstance(self.scene.attrs['orbit_number'], int))
        self.assertEqual(self.scene.attrs['orbit_number'], 12345)

    @mock.patch("level1c4pps.modis2pps_lib.Scene")
    def test_process_one_scene(self, mock_scene_class):
        """Test to set process_one_scene."""
        import level1c4pps.modis2pps_lib as modis2pps
        mock_scene_class.return_value = self.scene
        filename = modis2pps.process_one_scene("dummy", out_path='./level1c4pps/tests/', orbit_n='12345')
        self.assertEqual(os.path.basename(filename), "S_NWC_modis_eos2_12345_20090701T1201000Z_20090701T1201000Z.nc")
