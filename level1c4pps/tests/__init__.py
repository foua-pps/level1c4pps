import unittest

from level1c4pps.tests import (test_angles, test_seviri2pps, test_gac2pps)


def suite():
    """The global test suite.
    """
    mysuite = unittest.TestSuite()
    mysuite.addTests(test_angles.suite())
    mysuite.addTests(test_seviri2pps.suite())
    mysuite.addTests(test_gac2pps.suite())
    return mysuite


def load_tests(loader, tests, pattern):
    return suite()
