#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 level1c4pps developers
#
# This file is part of level1c4pps
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

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Nina Hakansson <nina.hakansson@smhi.se>
#   Adam.Dybbroe <adam.dybbroe@smhi.se>
#   Stephan Finkensieper <stephan.finkensieper@dwd.de>

# This program was developed by CMSAF to be used for the processing of
# CLAAS3.

"""Script to make seviri level1c in PPS-format with pytroll."""


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
    parser.add_argument('-o', '--out_dir', type=str, default='.',
                        required=False,
                        help="Output directory where to store level1c file.")
    parser.add_argument('--no-rotation', action='store_true',
                        help="Don't rotate images")
    parser.add_argument('-ne', '--nc_engine', type=str, nargs='?',
                        required=False, default='h5netcdf',
                        help="Engine for saving netcdf files netcdf4 or h5netcdf (default).")
    parser.add_argument('--use-nominal-time-in-filename', action='store_true',
                        help='Use nominal scan timestamps in output filename.')
    parser.add_argument('--no-sun-earth-distance-correction',
                        action='store_true',
                        help='Do not apply sun earth distance correction.')
    options = parser.parse_args()
    process_one_scan(
        options.files,
        out_path=options.out_dir,
        rotate=not options.no_rotation,
        engine=options.nc_engine,
        use_nominal_time_in_filename=options.use_nominal_time_in_filename,
        apply_sun_earth_distance_correction=not options.no_sun_earth_distance_correction
    )
