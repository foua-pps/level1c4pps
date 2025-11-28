#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2025 level1c4pps developers
#
# This file is part of level1c4pps.
#
# level1c4pps is free software: you can redistribute it and/or modify it
# under the terms of the Apache-2.0 license.
#
#

"""Unit tests for misc functions in __init__.py."""

import unittest
import xarray as xr
import level1c4pps
import numpy as np


class TestInit(unittest.TestCase):
    """Test functions in __init__.py."""

    def test_get_band_encoding(self):
        """Test get encoding."""
        ds = xr.DataArray([], attrs={'name': 'dummy'})
        self.assertRaises(ValueError, level1c4pps.get_band_encoding, ds,
                          None, None)

    def test_adjust_lons(self):
        """Test adjusted longitudes."""
        from level1c4pps import centered_modulus
        in_lons = xr.DataArray([340.0, 10.0, -22.0])
        out_lons = xr.DataArray([-20.0, 10.0, -22.0])
        in_lons_np = np.array([340.0, 10.0, -22.0])
        out_lons_np = np.array([-20.0, 10.0, -22.0])
        np.testing.assert_allclose(centered_modulus(in_lons),
                                   out_lons, rtol=0.00001)
        np.testing.assert_allclose(centered_modulus(in_lons_np),
                                   out_lons_np, rtol=0.00001)


def suite():
    """Create the test suite for test_init."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestInit))
    return mysuite
