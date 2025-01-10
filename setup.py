#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 level1c4pps developers

# Author(s): AUTHORS.md

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


install_requires = [
    "h5netcdf",
    "pyorbital",
    "pyspectral",
    "trollsift",
    "satpy[avhrr_l1b_eps,viirs_l1b,viirs_sdr,viirs_compact]!=0.38.*,!=0.39.*,!=0.40.*,!=0.41.*",
]
extras_requires = {
    "extra": "satpy[avhrr_l1b_gaclac,seviri_l1b_hrit,seviri_l1b_native,seviri_l1b_nc,vii_l1b_nc]"
}
NAME = "level1c4pps"
README = open('README.md', 'r').read()

setup(name=NAME,
      # version=version.__version__,
      description='Tools to convert various satellite level-1 data formats to NWCSAF/PPS level-1c format',
      long_description=README,
      author='NWCSAF PPS Team et al.',
      author_email='6160529+ninahakansson@users.noreply.github.com',
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
               'bin/mersi2pps.py',
               'bin/mersi22pps.py',
               'bin/viirs2pps.py',
               'bin/vgac2pps.py',
               'bin/slstr2pps.py',
               'bin/metimage2pps.py',
               'bin/eumgacfdr2pps.py',
               'bin/modis2pps.py',
               'bin/avhrr2pps.py'],
      data_files=[],
      zip_safe=False,
      use_scm_version=True,
      python_requires='>=3.7',
      install_requires=install_requires,
      extras_require=extras_requires,
      test_suite='level1c4pps.tests.suite',
      )
