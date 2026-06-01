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
import os
import unittest
from unittest import mock

import xarray as xr
from satpy import Scene

import level1c4pps.viirs2pps_lib as viirs2pps


class TestViirs2PPS(unittest.TestCase):
    """Test viirs2pps_lib."""

    def setUp(self):
        """Create a test scene."""
        self.scene = Scene()
        scene_dict = {}
        grid_data = [[1.0, 2.0], [3.0, 4.0]]
        all_keys = ['M05', 'M15', 'M12', 'I04'] + ['m_latitude', 'm_longitude'] + viirs2pps.ANGLE_NAMES
        for key in all_keys:
            scene_dict[key] = xr.DataArray(grid_data,
                                           dims=('y', 'x'),
                                           attrs={'name': key,
                                                  'id_tag': key})
        scene_dict['M05'].attrs = {'name': 'image0',
                                   'wavelength': [1, 2, 3, 'um'],
                                   'modifiers': ("sunz_correction", ),
                                   'id_tag': 'ch_r06'}
        scene_dict['M15'].attrs = {'name': 'image1',
                                   'id_tag': 'ch_tb11',
                                   'platform': "Suomi-NPP",
                                   'number_of_scans': 15,
                                   'rows_per_scan': 16,
                                   'wavelength': [1, 2, 3, 'um'],
                                   'start_time': dt.datetime(2009, 7, 1, 12, 1, 0),
                                   'end_time': dt.datetime(2009, 7, 1, 12, 1, 0),
                                   'orbit_number': 99999}
        scene_dict['M12'].attrs = {'name': 'image2',
                                   'id_tag': 'ch_tb37',
                                   'wavelength': [1, 2, 3, 'um']}
        scene_dict['I04'].attrs = {'name': 'image3',
                                   'id_tag': 'ch_tb11',
                                   'platform': "Suomi-NPP",
                                   'number_of_scans': 15,
                                   'rows_per_scan': 16,
                                   'wavelength': [1, 2, 3, 'um'],
                                   'start_time': dt.datetime(2009, 7, 1, 12, 1, 0),
                                   'end_time': dt.datetime(2009, 7, 1, 12, 1, 0),
                                   'orbit_number': 99999}
        for key in scene_dict:
            self.scene[key] = scene_dict[key]
        self.scene.load = mock.MagicMock()
        self.scene.resample = mock.MagicMock(return_value=self.scene)
        self.scene.attrs['sensor'] = ['viirs']

    def test_set_header_and_band_attrs(self):
        """Test to set header_and_band_attrs."""
        viirs2pps.set_header_and_band_attrs(self.scene)
        self.assertTrue(isinstance(self.scene.attrs['orbit_number'], int))
        self.assertTrue(self.scene["M05"].attrs['sun_zenith_angle_correction_applied'])

    @mock.patch("level1c4pps.viirs2pps_lib.check_file_exists")
    @mock.patch("level1c4pps.viirs2pps_lib.Scene")
    def test_process_one_scene(self, mock_scene_class, mock_check_file_exists):
        """Test to set process_one_scene."""
        mock_scene_class.return_value = self.scene
        filename = viirs2pps.process_one_scene("dummy", out_path='./level1c4pps/tests/', channel_selection="all")
        self.assertEqual(os.path.basename(filename), "S_NWC_viirs_npp_00000_20090701T1201000Z_20090701T1201000Z.nc")
        self.scene.load.assert_called_with(
            ['M01', 'M02', 'M03', 'M04', 'M05', 'M06', 'M07', 'M08',
             'M09', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16',
             'satellite_zenith_angle',
             'solar_zenith_angle',
             'satellite_azimuth_angle',
             'solar_azimuth_angle',
             'm_latitude',
             'm_longitude'], resolution=742)

    @mock.patch("level1c4pps.viirs2pps_lib.check_file_exists")
    @mock.patch("level1c4pps.viirs2pps_lib.Scene")
    def test_process_one_scene_iband(self, mock_scene_class, mock_check_file_exists):
        """Test to set process_one_scene."""
        mock_scene_class.return_value = self.scene
        filename = viirs2pps.process_one_scene("dummy", out_path='./level1c4pps/tests/', use_iband_res=True)
        self.assertEqual(os.path.basename(filename), "S_NWC_viirs_npp_00000_20090701T1201000Z_20090701T1201000Z.nc")
        self.scene.load.assert_called_with(['M09', 'M11', 'M14', 'M15', 'M16'], resolution=742)

    @mock.patch("level1c4pps.viirs2pps_lib.check_file_exists")
    @mock.patch("level1c4pps.viirs2pps_lib.Scene")
    def test_process_one_scene_compact(self, mock_scene_class, mock_check_file_exists):
        """Test to set process_one_scene."""
        mock_scene_class.return_value = self.scene
        filename = viirs2pps.process_one_scene("dummy", out_path='./level1c4pps/tests/', reader="viirs_compact")
        self.assertEqual(os.path.basename(filename), "S_NWC_viirs_npp_00000_20090701T1201000Z_20090701T1201000Z.nc")
        self.scene.load.assert_called_with(['M05', 'M07', 'M09', 'M10', 'M11'], modifiers=('sunz_corrected',))
