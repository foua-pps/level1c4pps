#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019, 2021 Pytroll

# Author(s):

#   Nina Håkansson <nina.hakansson@smhi.se>
#   Erik Johansson <erik.johansson@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Scripts to convert various satellite level-1 data formats to NWCSAF/PPS level-1c format."""

from setuptools import setup
from setuptools import find_packages

try:
    # HACK: https://github.com/pypa/setuptools_scm/issues/190#issuecomment-351181286
    # Stop setuptools_scm from including all repository files
    import setuptools_scm.integration
    setuptools_scm.integration.find_files = lambda _: []
except ImportError:
    pass

requires = ['satpy >= 0.19', 'pyorbital', 'trollsift', 'pyspectral', 'h5netcdf'],

NAME = "level1c4pps"
README = open('README.md', 'r').read()

setup(name=NAME,
      # version=version.__version__,
      description='Tools to convert various satellite level-1 data formats to NWCSAF/PPS level-1c format',
      long_description=README,
      author='Nina Hakansson',
      author_email='nina.hakansson@smhi.se',
      classifiers=["Development Status :: 3 - Alpha",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License v3 " +
                   "or later (GPLv3+)",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python",
                   "Topic :: Scientific/Engineering"],
      url="https://github.com/foua-pps/level1c4pps",
      packages=find_packages(),
      scripts=['bin/seviri2pps.py',
               'bin/gac2pps.py',
               'bin/mersi22pps.py',
               'bin/viirs2pps.py',
               'bin/slstr2pps.py',
               'bin/metimage2pps.py',
               'bin/eumgacfdr2pps.py',
               'bin/modis2pps.py',
               'bin/avhrr2pps.py'],
      data_files=[],
      zip_safe=False,
      use_scm_version=True,
      python_requires='>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*',
      install_requires=requires,
      test_suite='level1c4pps.tests.suite',
      )
