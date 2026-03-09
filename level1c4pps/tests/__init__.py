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

"""Test Initializer for level1c4pps."""

import unittest

from level1c4pps.tests import (test_angles, test_seviri2pps, test_gac2pps,
                               test_mersi2pps, test_modis2pps, test_slstr2pps,
                               test_viirs2pps, test_eumgacfdr2pps,
                               test_avhrr2pps, test_init)


def suite():
    """Test global test suite."""
    mysuite = unittest.TestSuite()
    mysuite.addTests(test_angles.suite())
    mysuite.addTests(test_seviri2pps.suite())
    mysuite.addTests(test_gac2pps.suite())
    mysuite.addTests(test_eumgacfdr2pps.suite())
    mysuite.addTests(test_mersi2pps.suite())
    mysuite.addTests(test_modis2pps.suite())
    mysuite.addTests(test_avhrr2pps.suite())
    mysuite.addTests(test_slstr2pps.suite())
    mysuite.addTests(test_viirs2pps.suite())
    mysuite.addTests(test_init.suite())
    return mysuite


def load_tests(loader, tests, pattern):
    """Load all tests."""
    return suite()
