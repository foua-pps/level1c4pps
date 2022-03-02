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

"""Unit tests for the eumgacfdr2pps_lib module."""

import netCDF4
import unittest
from datetime import datetime
try:
    from unittest import mock
except ImportError:
    import mock
from satpy import Scene

import level1c4pps.eumgacfdr2pps_lib as eumgacfdr2pps
import numpy as np


class TestEumgacfdr2PPS(unittest.TestCase):
    """Test eumgacfdr2pps_lib."""

    def setUp(self):
        """Create a test scene."""
        vis006 = mock.MagicMock(attrs={'name': 'image0',
                                       'wavelength': [1, 2, 3, 'um'],
                                       'id_tag': 'ch_r06'})
        ir_108 = mock.MagicMock(attrs={'name': 'image1',
                                       'id_tag': 'ch_tb11',
                                       'wavelength': [1, 2, 3, 'um'],
                                       'start_time': datetime.utcnow(),
                                       'end_time': datetime.utcnow(),
                                       'history': 'dummy',
                                       'platform_name': 'tirosn',
                                       'orbit_number': 99999})
        qual_f = mock.MagicMock(attrs={'name': 'qual_flags',
                                       'id_tag': 'qual_flags'})
        scan_t = mock.MagicMock(attrs={'name': 'scanline_timestamps'})
        self.scene = Scene()
        self.scene.attrs['sensor'] = ['avhrr-1', 'avhrr-2', 'avhrr-3']
        scene_dict = {'reflectance_channel_1': vis006, 'brightness_temperature_channel_4': ir_108,
                      'qual_flags': qual_f, 'acq_time': scan_t}
        for key in scene_dict:
            pps_name = scene_dict[key].attrs['name']
            self.scene[key] = scene_dict[key]
            self.scene[key].attrs['name'] = pps_name

    def test_get_encoding(self):
        """Test the encoding for EUMSAT GAC FDR."""
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
                                    'units': 'Milliseconds since 1970-01-01',
                                    'complevel': 4, '_FillValue': -1.0},
        }
        encoding = eumgacfdr2pps.get_encoding_gac(self.scene)
        self.assertDictEqual(encoding, encoding_exp)

    def test_set_header_and_band_attrs(self):
        """Test to set header_and_band_attrs."""
        eumgacfdr2pps.set_header_and_band_attrs(self.scene, orbit_n='12345')
        self.assertEqual(self.scene.attrs['orbit_number'], 12345)

    def test_process_one_file(self):
        """Test process one file for one example file."""
        # '1 11060U 78096A   80003.54792075  .00000937  00000-0  52481-3 0  2588\r\n',
        # '2 11060  98.9783 332.1605 0012789  88.8047 271.4583 14.11682873 63073\r\n']
        eumgacfdr2pps.process_one_file(
            './level1c4pps/tests/AVHRR-GAC_FDR_1C_N06_19810330T042358Z_19810330T060903Z_R_O_20200101T000000Z_0100.nc',
            out_path='./level1c4pps/tests/',
        )
        filename = './level1c4pps/tests/S_NWC_avhrr_noaa6_99999_19810330T0423582Z_19810330T0424032Z.nc'
        # written with hfnetcdf read with NETCDF4 ensure compatability
        pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')  # Check compatability implicitly
        for key in ['date_created', 'end_time', 'history', 'instrument',
                    'orbit_number', 'platform',
                    'sensor', 'source', 'start_time', 'Conventions']:
            if key not in pps_nc.__dict__.keys():
                print("Missing in attributes:", key)
            self.assertTrue(key in pps_nc.__dict__.keys())

        expected_vars = ['satzenith', 'azimuthdiff',
                         'satazimuth', 'sunazimuth', 'sunzenith',
                         'time', 'y', 'num_flags', 'lon', 'lat', 'qual_flags',
                         'image1', 'image3', 'image2', 'image5',
                         'midnight_line', 'overlap_free_end',
                         'overlap_free_start', 'x',
                         'equator_crossing_longitude',
                         'equator_crossing_time',
                         'scanline_timestamps', 'time_bnds']
        optional = ['bands_1d', 'acq_time']
        for var in optional:
            if var in pps_nc.variables.keys():
                expected_vars.append(var)
        self.assertEqual(sorted(pps_nc.variables.keys()),
                         sorted(expected_vars))

        np.testing.assert_almost_equal(pps_nc.variables['image1'].sun_earth_distance_correction_factor,
                                       0.9975245, decimal=4)


def suite():
    """Create the test suite for test_eumgacfdr2pps."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestEumgacfdr2PPS))

    return mysuite
