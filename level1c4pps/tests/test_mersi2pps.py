#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2025 level1c4pps developers
#
# This file is part of level1c4pps
#
# level1c4pps is free software: you can redistribute it and/or modify it
# under the terms of the Apache-2.0 license.
#
#

"""Unit tests for the merci2pps_lib module."""

import datetime as dt
import os
import unittest
from unittest import mock

import numpy as np
import xarray as xr
from satpy import Scene

import level1c4pps.mersi2pps_lib as mersi2pps


class TestMersi2PPS(unittest.TestCase):
    """Test mersi2pps_lib."""

    def setUp(self):
        """Create a test scene."""
        self.scene = Scene()
        scene_dict = {}
        grid_data = [[1.0, 2.0], [3.0, 4.0]]
        all_keys = ['3', '24'] + mersi2pps.GEOLOCATION_NAMES
        for key in all_keys:
            scene_dict[key] = xr.DataArray(grid_data,
                                           dims=('y', 'x'),
                                           attrs={'name': key,
                                                  'id_tag': key})
        scene_dict['3'].attrs = {'name': 'image0',
                                 'wavelength': [1, 2, 3, 'um'],
                                 'modifiers': ("sunz_correction", ),
                                 'id_tag': 'ch_r06'}
        scene_dict['24'].attrs = {'name': 'image1',
                                  'id_tag': 'ch_tb11',
                                  'platform': "FY-3F",
                                  'number_of_scans': 15,
                                  'rows_per_scan': 16,
                                  'wavelength': [1, 2, 3, 'um'],
                                  'start_time': dt.datetime(2009, 7, 1, 12, 1, 0),
                                  'end_time': dt.datetime(2009, 7, 1, 12, 1, 0),
                                  'orbit_number': 99999}
        for key in scene_dict:
            self.scene[key] = scene_dict[key]
        self.scene.load = mock.MagicMock
        self.scene.attrs['sensor'] = ['mersi3']

    def test_compose_filename(self):
        """Test compose filename for MERSI-2."""
        start_time = dt.datetime(2009, 7, 1, 12, 15)
        end_time = dt.datetime(2009, 7, 1, 12, 30)
        scene = mock.MagicMock(attrs={'start_time': start_time,
                                      'end_time': end_time,
                                      'orbit_number': '99999',
                                      'platform': 'Noaa19'})
        start_time = dt.datetime(2009, 7, 1, 12, 16)
        end_time = dt.datetime(2009, 7, 1, 12, 27)
        band = mock.MagicMock(attrs={'start_time': start_time,
                                     'end_time': end_time})
        fname_exp = '/out/path/S_NWC_mersi2_noaa19_99999_20090701T1216000Z_20090701T1227000Z.nc'
        fname = mersi2pps.compose_filename(scene, '/out/path', 'mersi2', band=band)
        self.assertEqual(fname, fname_exp)

    def test_set_header_and_band_attrs(self):
        """Test to set header_and_band_attrs."""
        mersi2pps.set_header_and_band_attrs(self.scene, self.scene['24'], orbit_n='12345')
        self.assertTrue(isinstance(self.scene.attrs['orbit_number'], int))
        self.assertEqual(self.scene.attrs['orbit_number'], 12345)

    def test_remove_broken_data(self):
        """Test remove broken data."""
        data = xr.Dataset(
            {
                '3': (('y', 'x'), [[100, 0, 100, 100]]),
                '20': (('y', 'x'), [[200, 0, 200, 200]]),
                '22': (('y', 'x'), [[300, 0, 300, 300]]),
            }
        )
        mersi2pps.remove_broken_data(data)
        expect = xr.Dataset(
            {
                '3': (('y', 'x'), [[100, 0, 100, 100]]),
                '20': (('y', 'x'), [[200, np.nan, 200, 200]]),
                '22': (('y', 'x'), [[300, np.nan, 300, 300]]),
            }
        )
        for band in expect:
            np.testing.assert_array_equal(data[band], expect[band])

    def test_get_sensor(self):
        """Test get sensor."""
        for filename, expect in [
            ('tf2019234102243.FY3D-X_MERSI_GEOQK_L1B.HDF', 'mersi-2'),
            ('tf2019234102243.FY3F-X_MERSI_GEOQK_L1B.HDF', 'mersi-3'),
            ('not_recognized_file', None),
        ]:
            sensor = mersi2pps.get_sensor(filename)
            self.assertEqual(sensor, expect)

    @mock.patch("level1c4pps.mersi2pps_lib.check_file_exists")
    @mock.patch("level1c4pps.mersi2pps_lib.Scene")
    def test_process_one_scene(self, mock_scene_class, mock_check_file_exists):
        """Test to set process_one_scene."""
        mock_scene_class.return_value = self.scene
        filename = mersi2pps.process_one_scene(
            ["tf2019234102243.FY3F-X_MERSI_GEOQK_L1B.HDF"], out_path='./level1c4pps/tests/', orbit_n='12345')
        self.assertEqual(os.path.basename(filename), "S_NWC_mersi3_fy3f_12345_20090701T1201000Z_20090701T1201000Z.nc")
