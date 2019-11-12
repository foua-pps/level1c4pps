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
# Author(s):

#   Stephan Finkensieper <stephan.finkensieper@dwd.de>
#   Nina Hakansson <nina.hakansson@smhi.se>

"""Unit tests for the gac2pps_lib module."""

import datetime as dt
import numpy as np
from pyresample.geometry import AreaDefinition
import unittest
try:
    from unittest import mock
except ImportError:
    import mock
import xarray as xr
from satpy import Scene

import level1c4pps.gac2pps_lib as gac2pps
import level1c4pps.calibration_coefs as calib


class TestGac2PPS(unittest.TestCase):
    """Test gac2pps_lib."""
    def setUp(self):
        """Create a test scene."""
        gac2pps.BANDNAMES = ['1', '4']
        vis006 = mock.MagicMock(attrs={'name': 'image0',
                                       'id_tag': 'ch_r06'})
        ir_108 = mock.MagicMock(attrs={'name': 'image1',
                                       'id_tag': 'ch_tb11',
                                       'platform_name': 'tirosn',
                                       'orbit_number': 99999})
        qual_f = mock.MagicMock(attrs={'name': 'qual_flags'})
        scan_t = mock.MagicMock(attrs={'name': 'scanline_timestamps'})
        self.scene = Scene()
        scene_dict = {'1': vis006, '4': ir_108, 'qual_flags': qual_f, 'scanline_timestamps': scan_t}
        for key in scene_dict:
            pps_name = scene_dict[key].attrs['name']
            self.scene[key] = scene_dict[key]
            self.scene[key].attrs['name'] = pps_name

    def test_get_encoding(self):
        enc_exp_angles = {'dtype': 'int16',
                          'scale_factor': 0.01,
                          'zlib': True,
                          'complevel': 4,
                          '_FillValue': -32767,
                          'add_offset': 0.0}
        enc_exp_coords = {'dtype': 'float32',
                          'zlib': True,
                          'complevel': 4,
                          '_FillValue': -999.0}
        encoding_exp = {
            'image0': {'dtype': 'int16',
                       'scale_factor': 0.01,
                       'zlib': True,
                       'complevel': 4,
                       '_FillValue': -32767,
                       'add_offset': 0.0},
            'image1': {'dtype': 'int16',
                       'scale_factor': 0.01,
                       '_FillValue': -32767,
                       'zlib': True,
                       'complevel': 4,
                       'add_offset': 273.15},
            'qual_flags':  {'dtype': 'int16', 'zlib': True,
                            'complevel': 4, '_FillValue': -32001.0},
            'scanline_timestamps': {'dtype': 'int64', 'zlib': True,
                                    'complevel': 4, '_FillValue': -1.0},
        }
        encoding = gac2pps.get_encoding_gac(self.scene)
        self.assertDictEqual(encoding, encoding_exp)

    def test_compose_filename(self):
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
        fname_exp = '/out/path/S_NWC_avhrr_noaa19_99999_20090701T1216000Z_20090701T1227000Z.nc'
        fname = gac2pps.compose_filename(scene, '/out/path', 'avhrr', band=band)
        self.assertEqual(fname, fname_exp)

    def test_set_header_and_band_attrs(self):
        gac2pps.set_header_and_band_attrs(self.scene)


def suite():
    """Create the test suite for test_gac2pps."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestGac2PPS))

    return mysuite
