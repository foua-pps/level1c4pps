#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 level1c4pps developers
#
# This file is part of level1c4pps
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
# Author(s):

#   Stephan Finkensieper <stephan.finkensieper@dwd.de>

"""Unit tests for the seviri2pps_lib module."""

import datetime as dt
import numpy as np
from pyresample.geometry import AreaDefinition
import unittest
from unittest import mock
import xarray as xr

import level1c4pps.seviri2pps_lib as seviri2pps
import level1c4pps.calibration_coefs as calib



class TestSeviri2PPS(unittest.TestCase):
    def test_rotate_band(self):
        """Test rotation of bands."""
        area = AreaDefinition(area_id='test',
                              description='test',
                              proj_id='test',
                              projection={'proj': 'geos', 'h': 12345},
                              width=3,
                              height=3,
                              area_extent=[1001, 1002, -1003, -1004])
        data = xr.DataArray(data=[[1, 2, 3],
                                  [4, 5, 6],
                                  [7, 8, 9]],
                            dims=('y', 'x'),
                            coords=[('y', [1.1, 0, -1.1]), ('x', [1, 0, -1])],
                            attrs={'area': area})
        scene = {'data': data}

        # Rotate
        seviri2pps.rotate_band(scene, 'data')

        # Check results
        self.assertTupleEqual(scene['data'].attrs['area'].area_extent,
                              (-1003, -1004, 1001, 1002))
        np.testing.assert_array_equal(scene['data']['x'], [-1, 0, 1])
        np.testing.assert_array_equal(scene['data']['y'], [-1.1, 0, 1.1])
        np.testing.assert_array_equal(scene['data'], [[9, 8, 7],
                                                      [6, 5, 4],
                                                      [3, 2, 1]])
        lons, lats = scene['data'].attrs['area'].get_lonlats()
        self.assertTrue(lons[0, 0] < 0)
        self.assertTrue(lons[0, 2] > 0)
        self.assertTrue(lats[0, 0] > 0)
        self.assertTrue(lons[2, 0] < 0)

    def test_get_lonlats(self):
        """Test getting lat/lon coordinates."""
        lons = np.array([1, 2, -1234, 1234], dtype=float)
        lats = np.array([-1234, 1234, 1, 2], dtype=float)
        area = mock.MagicMock()
        area.get_lonlats.return_value = lons, lats
        ds = mock.MagicMock(attrs={'area': area})

        lons_m, lats_m = seviri2pps.get_lonlats(ds)

        np.testing.assert_array_equal(lons_m, np.array([1, 2, np.nan, np.nan]))
        np.testing.assert_array_equal(lats_m, np.array([np.nan, np.nan, 1, 2]))

    @mock.patch('level1c4pps.seviri2pps_lib.get_mean_acq_time')
    @mock.patch('level1c4pps.seviri2pps_lib.sun_zenith_angle')
    @mock.patch('level1c4pps.seviri2pps_lib.get_alt_az')
    def test_get_solar_angles(self, get_alt_az, sun_zenith_angle,
                              get_mean_acq_time):
        """Test getting solar angles."""
        def sunz_patched(time, lon, lat):
            return time.astype(int) + lon + lat

        def alt_az_patched(time, lon, lat):
            return (time.astype(int) + lon + lat) * np.pi / 2

        get_alt_az.side_effect = alt_az_patched
        sun_zenith_angle.side_effect = sunz_patched
        get_mean_acq_time.return_value = xr.DataArray(np.array(
            ['1970-01-01 00:00:00.000000003',
             '1970-01-01 00:00:00.000000002',
             '1970-01-01 00:00:00.000000001',
             'NaT'], dtype='datetime64[ns]'))  # [3, 2, 1, -9E18] as int

        lons = np.array([[1, 2],
                         [3, 4],
                         [5, 6],
                         [0, 0]])
        lats = np.array([[-1, -2],
                         [-3, -4],
                         [-5, -6],
                         [0, 0]])
        suna_exp = np.array([[270, 270],
                             [180, 180],
                             [90, 90],
                             [np.nan, np.nan]])
        sunz_exp = np.array([[3, 3],
                             [2, 2],
                             [1, 1],
                             [np.nan, np.nan]])

        suna, sunz = seviri2pps.get_solar_angles('scene', lons=lons, lats=lats)
        np.testing.assert_array_equal(suna, suna_exp)
        np.testing.assert_array_equal(sunz, sunz_exp)

    @mock.patch('level1c4pps.seviri2pps_lib.get_observer_look')
    @mock.patch('level1c4pps.seviri2pps_lib.satpy.utils.get_satpos')
    def test_get_satellite_angles(self, get_satpos, get_observer_look):
        """Test getting satellite angles."""
        def get_observer_look_patched(lon, lat, alt, *args):
            if alt == 36000*1000:
                return None, 31  # > 30
            elif alt == 36000:
                return None, 22  # < 20
            else:
                return 'sata', 176

        get_observer_look.side_effect = get_observer_look_patched
        get_satpos.return_value = 'sat_lon', 'sat_lat', 12345678
        ds = mock.MagicMock(attrs={'start_time': 'start_time'})
        sata, satz = seviri2pps.get_satellite_angles(ds, 'lons', 'lats')
        self.assertEqual(sata, 'sata')
        self.assertEqual(satz, -86)
        get_observer_look.assert_called_with('sat_lon', 'sat_lat', 12345.678,
                                             'start_time', 'lons', 'lats', 0)

        # Height in km
        get_satpos.return_value = 'sat_lon', 'sat_lat', 36000
        self.assertRaises(seviri2pps.UnexpectedSatpyVersion,
                          seviri2pps.get_satellite_angles, ds, 'lons', 'lats')

        # pyorbital behaves unexpectedly
        get_satpos.return_value = 'sat_lon', 'sat_lat', 38001
        get_observer_look.reset_mock(side_effect=True)
        get_observer_look.return_value = None, 9999
        self.assertRaises(seviri2pps.UnexpectedSatpyVersion,
                          seviri2pps.get_satellite_angles, ds, 'lons', 'lats')

    def test_set_attrs(self):
        """Test setting scene attributes."""
        seviri2pps.BANDNAMES = ['VIS006', 'IR_108']
        vis006 = mock.MagicMock(attrs={})
        ir108 = mock.MagicMock(attrs={'platform_name': 'myplatform',
                                      'orbital_parameters': {'orb_a': 1,
                                                             'orb_b': 2},
                                      'georef_offset_corrected': True})
        scene_dict = {'VIS006': vis006, 'IR_108': ir108}
        scene = mock.MagicMock(attrs={})
        scene.__getitem__.side_effect = scene_dict.__getitem__

        seviri2pps.set_attrs(scene)
        self.assertEqual(scene['VIS006'].attrs['name'], 'image0')
        self.assertEqual(scene['VIS006'].attrs['id_tag'], 'ch_r06')
        self.assertEqual(scene['IR_108'].attrs['name'], 'image1')
        self.assertEqual(scene['IR_108'].attrs['id_tag'], 'ch_tb11')
        self.assertEqual(scene.attrs['orb_a'], 1)
        self.assertEqual(scene.attrs['orb_b'], 2)
        self.assertEqual(scene.attrs['georef_offset_corrected'], True)

    def test_get_mean_acq_time(self):
        """Test computation of mean scanline acquisition time."""
        seviri2pps.BANDNAMES = ['VIS006', 'IR_108']
        vis006 = xr.DataArray(
            data=[0, 0, 0],
            dims=('y', ),
            coords={'acq_time': ('y', [None,
                                       None,
                                       dt.datetime(2009, 7, 1, 12, 1, 0)])})
        ir_108 = xr.DataArray(
            data=[0, 0, 0],
            dims=('y', ),
            coords={'acq_time': ('y', [None,
                                       dt.datetime(2009, 7, 1, 12, 0, 30),
                                       dt.datetime(2009, 7, 1, 12, 1, 30)])})
        scene = {'VIS006': vis006, 'IR_108': ir_108}
        acq_exp = np.array(['NaT',
                            '2009-07-01 12:00:30',
                            '2009-07-01 12:01:15'], dtype='datetime64[s]')
        acq = seviri2pps.get_mean_acq_time(scene)
        np.testing.assert_array_equal(acq, acq_exp)

    @mock.patch('level1c4pps.seviri2pps_lib.get_mean_acq_time')
    def test_update_coords(self, get_mean_acq_time):
        """Test updating coordinates."""
        get_mean_acq_time.return_value = xr.DataArray([7, 8, 9], dims=('x',))
        seviri2pps.BANDNAMES = ['VIS006', 'IR_108']
        vis006 = xr.DataArray(data=[1, 2, 3],
                              dims=('x',),
                              coords={'acq_time': ('x', [0, 0, 0])},
                              attrs={'area': 'myarea',
                                     'start_time': dt.datetime(2009, 7, 1, 0)})
        ir_108 = xr.DataArray(data=[4, 5, 6],
                              dims=('x',),
                              coords={'acq_time': ('x', [0, 0, 0])},
                              attrs={'start_time': dt.datetime(2009, 7, 1, 1)})
        scene = {'VIS006': vis006.copy(), 'IR_108': ir_108.copy()}

        seviri2pps.update_coords(scene)

        for band in seviri2pps.BANDNAMES:
            self.assertEqual(scene[band].attrs['coordinates'], 'lon lat')
            np.testing.assert_array_equal(scene[band].coords['acq_time'].data,
                                          [7, 8, 9])

        np.testing.assert_array_equal(scene['VIS006'].data, vis006.data)
        np.testing.assert_array_equal(scene['IR_108'].data, ir_108.data)

        np.testing.assert_array_equal(scene['VIS006'].coords['time'].data,
                                      np.datetime64(dt.datetime(2009, 7, 1, 0)))
        np.testing.assert_array_equal(scene['IR_108'].coords['time'].data,
                                      np.datetime64(dt.datetime(2009, 7, 1, 1)))

    def test_add_ancillary_datasets(self):
        """Test adding ancillary datasets"""
        start_time = dt.datetime(2009, 7, 1, 0)
        end_time = dt.datetime(2009, 7, 1, 1)
        yvals = np.array([-1.0, 1.0])
        xvals = np.array([-1.1, 1.1])

        lons = np.array([[1.0, 2.0], [3.0, 4.0]])
        lats = np.array([[1.1, 2.1], [3.1, 4.1]])
        sunz = np.array([[1.2, 2.2], [3.2, 4.2]])
        satz = np.array([[1.3, 2.3], [3.3, 4.3]])
        azidiff = np.array([[1.4, 2.4], [3.4, 4.4]])

        ir_108 = xr.DataArray(data=np.array([[0.1, 0.2], [0.3, 0.4]]),
                              dims=('y', 'x'),
                              coords={'y': yvals,
                                      'x': xvals},
                              attrs={'start_time': start_time,
                                     'end_time': end_time,
                                     'orbital_parameters': 'orb_params',
                                     'georef_offset_corrected': True})
        scene = {'IR_108': ir_108}
        seviri2pps.add_ancillary_datasets(scene, lons=lons, lats=lats,
                                          sunz=sunz, satz=satz,
                                          azidiff=azidiff)

        # Test lon/lat
        np.testing.assert_array_equal(scene['lon'].data, lons)
        self.assertEqual(scene['lon'].attrs['units'], 'degrees_east')

        np.testing.assert_array_equal(scene['lat'].data, lats)
        self.assertEqual(scene['lat'].attrs['units'], 'degrees_north')

        # Test angles
        np.testing.assert_array_equal(scene['sunzenith'].data, sunz)
        self.assertEqual(scene['sunzenith'].attrs['name'], 'image11')

        np.testing.assert_array_equal(scene['satzenith'].data, satz)
        self.assertEqual(scene['satzenith'].attrs['name'], 'image12')

        np.testing.assert_array_equal(scene['azimuthdiff'].data, azidiff)
        self.assertEqual(scene['azimuthdiff'].attrs['name'], 'image13')

        for angle in ['azimuthdiff', 'satzenith', 'sunzenith']:
            self.assertTupleEqual(scene[angle].dims, ('y', 'x'))
            np.testing.assert_array_equal(scene[angle].coords['x'].data, xvals)
            np.testing.assert_array_equal(scene[angle].coords['y'].data, yvals)
            np.testing.assert_array_equal(scene[angle].coords['time'].data,
                                          np.datetime64(start_time))
            self.assertEqual(scene[angle].attrs['units'], 'degree')

        # Test common properties
        for name in ['lon', 'lat', 'azimuthdiff', 'satzenith', 'sunzenith']:
            self.assertTupleEqual(scene[name].dims, ('y', 'x'))
            np.testing.assert_array_equal(scene[name].coords['x'].data, xvals)
            np.testing.assert_array_equal(scene[name].coords['y'].data, yvals)
            self.assertEqual(scene[name].attrs['start_time'], start_time)
            self.assertEqual(scene[name].attrs['end_time'], end_time)

    def test_compose_filename(self):
        start_time = dt.datetime(2009, 7, 1, 12, 15)
        end_time = dt.datetime(2009, 7, 1, 12, 30)
        scene = mock.MagicMock(attrs={'start_time': start_time,
                                      'end_time': end_time,
                                      'platform': 'Meteosat-9'})
        fname_exp = '/out/path/S_NWC_seviri_meteosat9_99999_20090701T1215000Z_20090701T1230000Z.nc'
        fname = seviri2pps.compose_filename(scene, '/out/path')
        self.assertEqual(fname, fname_exp)

    def test_get_encoding(self):
        seviri2pps.BANDNAMES = ['VIS006', 'IR_108']
        vis006 = mock.MagicMock(attrs={'name': 'image0'})
        ir_108 = mock.MagicMock(attrs={'name': 'image1'})
        scene = {'VIS006': vis006, 'IR_108': ir_108}
        enc_exp_angles = {'dtype': 'int16',
                          'scale_factor': 0.01,
                          'zlib': True,
                          'complevel': 4,
                          '_FillValue': -32767,
                          'add_offset': 0.0,
                          'chunksizes': (1, 512, 3712)}
        enc_exp_coords = {'dtype': 'float32',
                          'zlib': True,
                          'complevel': 4,
                          '_FillValue': -999.0,
                          'chunksizes': (512, 3712)}
        encoding_exp = {
            'image0': {'dtype': 'int16',
                       'scale_factor': 0.01,
                       'zlib': True,
                       'complevel': 4,
                       '_FillValue': -32767,
                       'add_offset': 0.0,
                       'chunksizes': (1, 512, 3712)},
            'image1': {'dtype': 'int16',
                        'scale_factor': 0.01,
                        '_FillValue': -32767,
                        'zlib': True,
                        'complevel': 4,
                        'add_offset': 273.15,
                        'chunksizes': (1, 512, 3712)},
            'image11': enc_exp_angles,
            'image12': enc_exp_angles,
            'image13': enc_exp_angles,
            'lon': enc_exp_coords,
            'lat': enc_exp_coords
        }
        encoding = seviri2pps.get_encoding(scene)
        self.assertDictEqual(encoding, encoding_exp)

    def test_get_header_attrs(self):
        start_time = dt.datetime(2009, 7, 1, 12, 15)
        end_time = dt.datetime(2009, 7, 1, 12, 30)
        scene = mock.MagicMock(attrs={'foo': 'bar',
                                      'start_time': start_time,
                                      'end_time': end_time,
                                      'sensor': 'SEVIRI'})
        header_attrs_exp = {
            'foo': 'bar',
            'start_time': start_time,
            'end_time': end_time,
        }
        header_attrs = seviri2pps.get_header_attrs(scene)
        self.assertDictEqual(header_attrs, header_attrs_exp)


class TestCalibration(unittest.TestCase):
    def test_get_calibration_for_date(self):
        """Test MODIS-intercalibrated gain and offset for specific date."""
        coefs = calib.get_calibration_for_date(
            platform='MSG3', date=dt.date(2018, 1, 18))
        REF = {
            'VIS006': {'gain': 0.023689275200000002, 'offset': -1.2081530352},
            'VIS008': {'gain': 0.029757990399999996,
                       'offset': -1.5176575103999999},
            'IR_016': {'gain': 0.0228774688, 'offset': -1.1667509087999999}}
        for channel in REF.keys():
            self.assertEqual(coefs[channel]['gain'], REF[channel]['gain'])
            self.assertEqual(coefs[channel]['offset'], REF[channel]['offset'])

    def test_get_calibration_for_time(self):
        """Test MODIS-intercalibrated gain and offset for specific time."""
        REF = {
            ('MSG1', dt.datetime(2005, 1, 18, 0, 0)): {
                'VIS006': {'gain': 0.0250354716,
                           'offset': -1.2768090516000001},
                'VIS008': {'gain': 0.0315626684,
                           'offset': -1.6096960884},
                'IR_016': {'gain': 0.022880986,
                           'offset': -1.166930286}},
            ('MSG2', dt.datetime(2010, 1, 18, 0, 0)): {
                'VIS006': {'gain': 0.021964051999999998,
                           'offset': -1.120166652},
                'VIS008': {'gain': 0.027548445,
                           'offset': -1.404970695},
                'IR_016': {'gain': 0.021576766,
                           'offset': -1.100415066}},
            ('MSG3', dt.datetime(2018, 1, 18, 0, 0)): {
                 'VIS006': {'gain': 0.023689275200000002,
                            'offset': -1.2081530352},
                 'VIS008': {'gain': 0.029757990399999996,
                            'offset': -1.5176575103999999},
                 'IR_016': {'gain': 0.0228774688,
                            'offset': -1.1667509087999999}},
            ('MSG4', dt.datetime(2019, 1, 18, 0, 0)): {
                'VIS006': {'gain': 0.0230358454,
                           'offset': -1.1748281154},
                'VIS008': {'gain': 0.0292054763,
                           'offset': -1.4894792913000001},
                'IR_016': {'gain': 0.022189932800000003,
                           'offset': -1.1316865728}}
        }
        for (platform, time), ref in REF.items():
            coefs = calib.get_calibration_for_time(platform=platform,
                                                   time=time)
            for channel in ref.keys():
                np.testing.assert_allclose(coefs[channel]['gain'],
                                           ref[channel]['gain'])
                np.testing.assert_allclose(coefs[channel]['offset'],
                                           ref[channel]['offset'])

    def test_get_calibration(self):
        """Test MODIS-intercalibrated for date and time."""
        coefs1 = calib.get_calibration_for_time(
            platform='MSG3', time=dt.datetime(2018, 1, 18, 23, 59))
        coefs2 = calib.get_calibration_for_date(
            platform='MSG3', date=dt.date(2018, 1, 19))
        for channel in coefs1.keys():
            self.assertAlmostEqual(coefs1[channel]['gain'],
                                   coefs2[channel]['gain'],
                                   delta=10e-8)
            self.assertAlmostEqual(coefs1[channel]['offset'],
                                   coefs2[channel]['offset'],
                                   delta=10e-8)


def suite():
    """Create the test suite for test_seviri2pps."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestSeviri2PPS))
    mysuite.addTest(loader.loadTestsFromTestCase(TestCalibration))

    return mysuite