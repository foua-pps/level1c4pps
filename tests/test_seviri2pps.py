"""Unit tests for the seviri2pps_lib module."""

import unittest
import numpy as np
import xarray as xr


import level1c4pps.seviri2pps_lib as seviri2pps
from pyresample.geometry import AreaDefinition


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
