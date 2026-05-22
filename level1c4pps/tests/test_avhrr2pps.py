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

"""Unit tests for the avhrr2pps_lib module."""

import datetime as dt
import unittest
import os
from unittest import mock
from satpy import Scene
import xarray as xr

import level1c4pps.avhrr2pps_lib as avhrr2pps


class TestAvhrr2PPS(unittest.TestCase):
    """Test avhrr2pps_lib."""

    def setUp(self):
        """Create a test scene."""
        self.scene = Scene()
        scene_dict = {}
        for key in avhrr2pps.GEOLOCATION_NAMES_AAPP + ['1', '4']:
            scene_dict[key] = xr.DataArray([[1.0, 2.0],
                                            [3.0, 4.0]],
                                           dims=('y', 'x'),
                                           attrs={'name': key,
                                                  'id_tag': key})
        scene_dict['1'].attrs = {'name': 'image0',
                                  'wavelength': [1, 2, 3, 'um'],
                                  'id_tag': 'ch_r06'}
        scene_dict['4'].attrs = {'name': 'image1',
                                  'id_tag': 'ch_tb11',
                                  'platform': "NOAA-19",
                                  'wavelength': [1, 2, 3, 'um'],
                                  'start_time': dt.datetime(2009, 7, 1, 12, 1, 0),
                                  'end_time': dt.datetime(2009, 7, 1, 12, 1, 0),
                                  'platform_name': '',
                                  'orbit_number': 99999}
        for key in scene_dict:
            self.scene[key] = scene_dict[key]        
        self.scene.attrs['sensor'] = ['avhrr']



    def test_compose_filename(self):
        """Test compose filename for AVHRR."""
        start_time = dt.datetime(2009, 7, 1, 12, 15)
        end_time = dt.datetime(2009, 7, 1, 12, 30)
        scene = mock.MagicMock(attrs={'start_time': start_time,
                                      'end_time': end_time,
                                      'orbit_number': '99999',
                                      'platform': 'Metop-B'})
        start_time = dt.datetime(2009, 7, 1, 12, 16)
        end_time = dt.datetime(2009, 7, 1, 12, 27)
        band = mock.MagicMock(attrs={'start_time': start_time,
                                     'end_time': end_time})
        fname_exp = '/out/path/S_NWC_avhrr_metopb_99999_20090701T1216000Z_20090701T1227000Z.nc'
        fname = avhrr2pps.compose_filename(scene, '/out/path', 'avhrr', band=band)
        self.assertEqual(fname, fname_exp)

    def test_set_header_and_band_attrs(self):
        """Test to set header_and_band_attrs."""
        avhrr2pps.set_header_and_band_attrs(self.scene, orbit_n='12345')
        self.assertTrue(isinstance(self.scene.attrs['orbit_number'], int))
        self.assertEqual(self.scene.attrs['orbit_number'], 12345)
        
    @mock.patch("level1c4pps.avhrr2pps_lib.load_data")
    def test_process_one_scene(self, mock_load):
        """Test to set process_one_scene."""
        import level1c4pps.avhrr2pps_lib as avhrr2pps
        mock_load.return_value = self.scene
        filename = avhrr2pps.process_one_scene("dummy", out_path='./level1c4pps/tests/', orbit_n='12345')
        self.assertEqual(os.path.basename(filename), "S_NWC_avhrr_noaa19_12345_20090701T1201000Z_20090701T1201000Z.nc")
