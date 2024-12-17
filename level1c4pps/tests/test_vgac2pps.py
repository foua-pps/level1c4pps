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

"""Unit tests for the vgac2pps_lib module."""

import netCDF4
import unittest
from datetime import datetime, timezone
try:
    from unittest import mock
except ImportError:
    import mock
from satpy import Scene

import level1c4pps.vgac2pps_lib as vgac2pps
import numpy as np
import pytest

try:
    import sbafs_ann
    no_sbaf_module = False
except ModuleNotFoundError:
    no_sbaf_module = True


class TestVgac2PPS(unittest.TestCase):
    """Test vgac2pps_lib."""

    def setUp(self):
        """Create a test scene."""
        vis006 = mock.MagicMock(attrs={'name': 'image0',
                                       'wavelength': [1, 2, 3, 'um'],
                                       'id_tag': 'ch_r06'})
        ir_108 = mock.MagicMock(attrs={'name': 'image1',
                                       'id_tag': 'ch_tb11',
                                       'wavelength': [1, 2, 3, 'um'],
                                       'start_time': datetime.now(timezone.utc),
                                       'end_time': datetime.now(timezone.utc),
                                       'history': 'dummy',
                                       'platform_name': 'tirosn',
                                       'orbit_number': 99999})
        qual_f = mock.MagicMock(attrs={'name': 'qual_flags',
                                       'id_tag': 'qual_flags'})
        scan_t = mock.MagicMock(attrs={'name': 'scanline_timestamps'})
        self.scene = Scene()
        self.scene.attrs['sensor'] = ['avhrr-1', 'avhrr-2', 'avhrr-3']
        scene_dict = {'M05': vis006, 'M15': ir_108,
                      'qual_flags': qual_f, 'acq_time': scan_t}
        for key in scene_dict:
            pps_name = scene_dict[key].attrs['name']
            self.scene[key] = scene_dict[key]
            self.scene[key].attrs['name'] = pps_name

    def test_get_encoding(self):
        """Test the encoding for VGAC."""
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
            'scanline_timestamps': {'dtype': 'int64',
                                    'zlib': True,
                                    'units': 'milliseconds since 1970-01-01',
                                    'complevel': 4,
                                    '_FillValue': -1.0},
        }
        encoding = vgac2pps.get_encoding_viirs(self.scene)
        self.assertDictEqual(encoding, encoding_exp)

    def test_set_header_and_band_attrs(self):
        """Test to set header_and_band_attrs."""
        vgac2pps.set_header_and_band_attrs(self.scene, orbit_n='12345')
        self.assertEqual(self.scene.attrs['orbit_number'], 12345)

    def test_process_one_scene(self):
        """Test process one scene for one example file."""
        vgac2pps.process_one_scene(
            ['./level1c4pps/tests/VGAC_VJ102MOD_A2018305_1042_n004946_K005.nc'],
            out_path='./level1c4pps/tests/'
        )
        filename = './level1c4pps/tests/S_NWC_viirs_noaa20_00000_20181101T1042080Z_20181101T1224090Z.nc'
        # written with hfnetcdf read with NETCDF4 ensure compatability
        pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')  # Check compatability implicitly

        for key in ['start_time', 'end_time', 'history', 'instrument',
                    'orbit_number', 'platform',
                    'sensor', 'source']:
            if key not in pps_nc.__dict__.keys():
                print("Missing in attributes:", key)
            self.assertTrue(key in pps_nc.__dict__.keys())

        expected_vars = ['satzenith', 'azimuthdiff',
                         'satazimuth', 'sunazimuth', 'sunzenith',
                         'lon', 'lat',
                         'image1', 'image2', 'image3', 'image4', 'image5',
                         'image6', 'image7', 'image8', 'image9',
                         'scanline_timestamps', 'time', 'time_bnds']
        for var in expected_vars:
            self.assertTrue(var in pps_nc.variables.keys())

        np.testing.assert_almost_equal(pps_nc.variables['image1'].sun_earth_distance_correction_factor,
                                       1.0, decimal=4)
        self.assertTrue(pps_nc.variables['image1'].sun_earth_distance_correction_applied)
        assert (
            pps_nc.variables["scanline_timestamps"].units == "milliseconds since 1970-01-01"
            and pps_nc.variables["scanline_timestamps"].dtype == "int"
        )

    def test_process_one_scene_n19(self):
        """Test process one scene for one example file."""

        vgac2pps.process_one_scene(
            ['./level1c4pps/tests/VGAC_VJ102MOD_A2018305_1042_n004946_K005.nc'],
            out_path='./level1c4pps/tests/',
            noaa19_sbaf_version='v6'
        )
        filename = './level1c4pps/tests/S_NWC_avhrr_vgac20_00000_20181101T1042080Z_20181101T1224090Z.nc'
        filename_viirs = './level1c4pps/tests/S_NWC_viirs_noaa20_00000_20181101T1042080Z_20181101T1224090Z.nc'
        # written with hfnetcdf read with NETCDF4 ensure compatability
        pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')  # Check compatability implicitly
        pps_nc_viirs = netCDF4.Dataset(filename_viirs, 'r', format='NETCDF4')  # Check compatability implicitly

        for key in ['start_time', 'end_time', 'history', 'instrument',
                    'orbit_number', 'platform',
                    'sensor', 'source']:
            if key not in pps_nc.__dict__.keys():
                print("Missing in attributes:", key)
            self.assertTrue(key in pps_nc.__dict__.keys())

        expected_vars = ['satzenith', 'azimuthdiff',
                         'satazimuth', 'sunazimuth', 'sunzenith',
                         'lon', 'lat',
                         'image1', 'image2', 'image3', 'image4', 'image5',
                         'scanline_timestamps', 'time', 'time_bnds']

        for var in expected_vars:
            self.assertTrue(var in pps_nc.variables.keys())

        np.testing.assert_almost_equal(pps_nc.variables['image1'].sun_earth_distance_correction_factor,
                                       1.0, decimal=4)

        np.testing.assert_equal(pps_nc.__dict__["platform"], "vgac20")
        self.assertTrue(np.abs(pps_nc.variables['image1'][0, 0, 0] - pps_nc_viirs.variables['image1'][0, 0, 0]) > 0.01)

    @unittest.skipIf(no_sbaf_module, "Install sbafs_ann to test NN-SBAFS.")
    def test_process_one_scene_n19_nn(self):
        """Test process one scene for one example file."""
        import sbafs_ann

        vgac2pps.process_one_scene(
            ['./level1c4pps/tests/VGAC_VJ102MOD_A2018305_1042_n004946_K005.nc'],
            out_path='./level1c4pps/tests/',
            noaa19_sbaf_version='NN_v4'
        )
        filename = './level1c4pps/tests/S_NWC_avhrr_vgac20_00000_20181101T1042080Z_20181101T1224090Z.nc'
        filename_viirs = './level1c4pps/tests/S_NWC_viirs_noaa20_00000_20181101T1042080Z_20181101T1224090Z.nc'
        # written with hfnetcdf read with NETCDF4 ensure compatability
        pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')  # Check compatability implicitly
        pps_nc_viirs = netCDF4.Dataset(filename_viirs, 'r', format='NETCDF4')  # Check compatability implicitly

        for key in ['start_time', 'end_time', 'history', 'instrument',
                    'orbit_number', 'platform',
                    'sensor', 'source']:
            if key not in pps_nc.__dict__.keys():
                print("Missing in attributes:", key)
            self.assertTrue(key in pps_nc.__dict__.keys())

        expected_vars = ['satzenith', 'azimuthdiff',
                         'satazimuth', 'sunazimuth', 'sunzenith',
                         'lon', 'lat',
                         'image1', 'image2', 'image3', 'image4', 'image5',
                         'scanline_timestamps', 'time', 'time_bnds']

        for var in expected_vars:
            self.assertTrue(var in pps_nc.variables.keys())

        np.testing.assert_almost_equal(pps_nc.variables['image1'].sun_earth_distance_correction_factor,
                                       1.0, decimal=4)

        np.testing.assert_equal(pps_nc.__dict__["platform"], "vgac20")
        self.assertTrue(np.abs(pps_nc.variables['image1'][0, 0, 400] -
                        pps_nc_viirs.variables['image1'][0, 0, 400]) > 0.01)

    def test_process_one_scene_midnight(self):
        """Test process one scene for one example file."""

        vgac2pps.process_one_scene(
            ['./level1c4pps/tests/VGAC_VNPP02MOD_A2012365_2304_n06095_K005.nc'],
            out_path='./level1c4pps/tests/'
        )
        filename = './level1c4pps/tests/S_NWC_viirs_npp_00000_20121230T2359563Z_20121230T2359599Z.nc'
        # written with hfnetcdf read with NETCDF4 ensure compatability
        pps_nc = netCDF4.Dataset(filename, 'r', format='NETCDF4')  # Check compatability implicitly

        for key in ['start_time', 'end_time', 'history', 'instrument',
                    'orbit_number', 'platform',
                    'sensor', 'source']:
            if key not in pps_nc.__dict__.keys():
                print("Missing in attributes:", key)
            self.assertTrue(key in pps_nc.__dict__.keys())

        expected_vars = ['satzenith', 'azimuthdiff',
                         'satazimuth', 'sunazimuth', 'sunzenith',
                         'lon', 'lat',
                         'image1', 'image2', 'image3', 'image4', 'image5',
                         'image6', 'image7', 'image8', 'image9',
                         'scanline_timestamps', 'time', 'time_bnds']
        for var in expected_vars:
            self.assertTrue(var in pps_nc.variables.keys())

        print(pps_nc.variables['image1'].shape)

        np.testing.assert_equal(pps_nc.variables['image1'].shape, (1, 7, 801))
