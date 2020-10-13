#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 level1c4pps developers
#
# This file is part of level1c4pps.
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

#   Nina Hakansson <nina.hakansson@smhi.se>

"""Test Initializer for level1c4pps."""

import unittest

from level1c4pps.tests import (test_angles, test_seviri2pps, test_gac2pps,
                               test_mersi22pps,  test_modis2pps, test_slstr2pps,
                               test_viirs2pps, test_eumgacfdr2pps,
                               test_avhrr2pps, test_init)


def suite():
    """Test global test suite."""
    mysuite = unittest.TestSuite()
    mysuite.addTests(test_angles.suite())
    mysuite.addTests(test_seviri2pps.suite())
    mysuite.addTests(test_gac2pps.suite())
    mysuite.addTests(test_eumgacfdr2pps.suite())
    mysuite.addTests(test_mersi22pps.suite())
    mysuite.addTests(test_modis2pps.suite())
    mysuite.addTests(test_avhrr2pps.suite())
    mysuite.addTests(test_slstr2pps.suite())
    mysuite.addTests(test_viirs2pps.suite())
    mysuite.addTests(test_init.suite())
    return mysuite


def load_tests(loader, tests, pattern):
    """Load all tests."""
    return suite()
