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
import netCDF4
import unittest
from datetime import datetime
try:
    from unittest import mock
except ImportError:
    import mock
from satpy import Scene

import level1c4pps.gac2pps_lib as gac2pps
import numpy as np


class TestGac2PPS(unittest.TestCase):
    """Test gac2pps_lib."""

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
                                       'platform_name': 'tirosn',
                                       'orbit_number': 99999})
        qual_f = mock.MagicMock(attrs={'name': 'qual_flags',
                                       'id_tag': 'qual_flags'})
        scan_t = mock.MagicMock(attrs={'name': 'scanline_timestamps'})
        self.scene = Scene()
        self.scene.attrs['sensor'] = ['avhrr-1', 'avhrr-2', 'avhrr-3']
        scene_dict = {'1': vis006, '4': ir_108, 'qual_flags': qual_f, 'scanline_timestamps': scan_t}
        for key in scene_dict:
            pps_name = scene_dict[key].attrs['name']
            self.scene[key] = scene_dict[key]
            self.scene[key].attrs['name'] = pps_name

    def test_get_encoding(self):
        """Test the encoding for GAC."""
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
        encoding = gac2pps.get_encoding_gac(self.scene)
        self.assertDictEqual(encoding, encoding_exp)

    def test_compose_filename(self):
        """Test compose filename for GAC."""
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
        """Test to set header_and_band_attrs."""
        gac2pps.set_header_and_band_attrs(self.scene, orbit_n='12345')
        self.assertEqual(self.scene.attrs['orbit_number'], 12345)

    def test_process_one_file(self):
        """Test process one file for one example file."""
        # '1 11060U 78096A   80003.54792075  .00000937  00000-0  52481-3 0  2588\r\n',
        # '2 11060  98.9783 332.1605 0012789  88.8047 271.4583 14.11682873 63073\r\n']
        tle_dir = './level1c4pps/tests/'
        tle_name = 'TLE_tirosn.txt'
        gac2pps.process_one_file(
            './level1c4pps/tests/NSS.GHRR.TN.D80003.S1147.E1332.B0630506.GC',
            out_path='./level1c4pps/tests/',
            reader_kwargs={
                'tle_dir': tle_dir,
                'tle_name': tle_name
            })
        filename = './level1c4pps/tests/S_NWC_avhrr_tirosn_99999_19800103T1147154Z_19800103T1147229Z.nc'
        pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')
        for key in ['date_created', 'end_time', 'history', 'instrument',
                    'orbit_number', 'platform', 'platform_name',
                    'sensor', 'source', 'start_time', 'Conventions',
                    'version_level1c4pps', 'version_level1c4pps_satpy']:
            self.assertTrue(key in sorted(pps_nc.__dict__.keys()))
        expected_vars = ['satzenith', 'azimuthdiff', 'satazimuth', 'sunazimuth', 'sunzenith',
                         'time', 'y', 'num_flags', 'lon', 'lat', 'qual_flags',
                         'image1', 'image5', 'image3', 'image2',
                         'scanline_timestamps', 'time_bnds']
        optional = ['bands_1d', 'acq_time']
        for var in optional:
            if var in pps_nc.variables.keys():
                expected_vars.append(var)
        self.assertEqual(sorted(pps_nc.variables.keys()),
                         sorted(expected_vars))

        np.testing.assert_almost_equal(pps_nc.variables['image1'].sun_earth_distance_correction_factor,
                                       0.9666, decimal=4)

    def test_process_one_file_netcdf4(self):
        """Test process one file for one example file with the netcdf4 engine."""
        # '1 11060U 78096A   80003.54792075  .00000937  00000-0  52481-3 0  2588\r\n',
        # '2 11060  98.9783 332.1605 0012789  88.8047 271.4583 14.11682873 63073\r\n']
        tle_dir = './level1c4pps/tests/'
        tle_name = 'TLE_tirosn.txt'
        gac2pps.process_one_file(
            './level1c4pps/tests/NSS.GHRR.TN.D80003.S1147.E1332.B0630506.GC',
            out_path='./level1c4pps/tests/',
            reader_kwargs={
                'tle_dir': tle_dir,
                'tle_name': tle_name
            },
            engine='netcdf4')
        filename = './level1c4pps/tests/S_NWC_avhrr_tirosn_99999_19800103T1147154Z_19800103T1147229Z.nc'
        pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')

        expected_vars = ['satzenith', 'azimuthdiff', 'satazimuth', 'sunazimuth', 'sunzenith',
                         'time', 'y', 'num_flags', 'lon', 'lat', 'qual_flags',
                         'image1', 'image5', 'image3', 'image2',
                         'scanline_timestamps', 'time_bnds']
        optional = ['bands_1d', 'acq_time']
        for var in optional:
            if var in pps_nc.variables.keys():
                expected_vars.append(var)
        self.assertEqual(sorted(pps_nc.variables.keys()),
                         sorted(expected_vars))

        np.testing.assert_almost_equal(pps_nc.variables['image1'].sun_earth_distance_correction_factor,
                                       0.9666, decimal=4)


def suite():
    """Create the test suite for test_gac2pps."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestGac2PPS))

    return mysuite
