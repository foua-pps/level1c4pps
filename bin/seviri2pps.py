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

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Nina Hakansson <nina.hakansson@smhi.se>
#   Adam.Dybbroe <adam.dybbroe@smhi.se>

# This program was developed by CMSAF to be used for the processing of
# CLAAS3.

"""Script to make seviri level1c in PPS-format with pytroll"""


import argparse
from level1c4pps.seviri2pps_lib import process_one_scan

# -----------------------------------------------------------------------------
# Main:
if __name__ == "__main__":
    """ Create PPS-format level1c data 
    From a list of hirt files hrit create a level1c file for pps.
    """
    # python3 seviri2pps.py file1 file2 ... fileN -o output
    parser = argparse.ArgumentParser(
        description=('Script to produce a PPS-level1c file for a list of '
                     'SEVIRI hrit files.'))
    parser.add_argument('files', metavar='fileN', type=str, nargs='+',
                        help='List of hrit files to process for one scan')
    parser.add_argument('-o', '--out_dir', type=str, nargs='?',
                        required=False,
                        help="Output directory where to store level1c file.")
    parser.add_argument('-b', '--buggy_satz', const=True, nargs='?',
                        required=False,
                        help="Create buggy satellite zenith angle data")
    options = parser.parse_args()
    process_buggy_satellite_zenith_angles = False
    if options.buggy_satz:
        process_buggy_satellite_zenith_angles = True
    process_one_scan(options.files, options.out_dir,
                     process_buggy_satellite_zenith_angles)
