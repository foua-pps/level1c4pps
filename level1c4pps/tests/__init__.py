import unittest

from level1c4pps.tests import (test_angles, test_seviri2pps)


def suite():
    """The global test suite.
    """
    mysuite = unittest.TestSuite()
    mysuite.addTests(test_angles.suite())
    mysuite.addTests(test_seviri2pps.suite())
    return mysuite


def load_tests(loader, tests, pattern):
    return suite()
