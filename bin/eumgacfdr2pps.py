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

"""Script to make EUMETSAT GAC level1c in PPS-format with pytroll."""

import argparse
from level1c4pps.eumgacfdr2pps_lib import process_one_file


# -----------------------------------------------------------------------------
# Main:
if __name__ == "__main__":
    """ Create PPS-format level1c data
    From a list of hirt files hrit create a level1c file for pps.
    """
    parser = argparse.ArgumentParser(
        description=('Script to produce a PPS-level1c file for a '
                     'GAC.'))
    parser.add_argument('file', type=str,
                        help='GAC file to process')
    parser.add_argument('-o', '--out_dir', type=str, nargs='?',
                        required=False, default='.',
                        help="Output directory where to store level1c file.")
    parser.add_argument('-sl', '--start_line', type=int, nargs='?',
                        required=False, default=None,
                        help="Use only part of data, start with start_line.")
    parser.add_argument('-el', '--end_line', type=int, nargs='?',
                        required=False, default=None,
                        help="Use only part of data, start with end_line.")
    parser.add_argument('-ne', '--nc_engine', type=str, nargs='?',
                        required=False, default='h5netcdf',
                        help="Engine for saving netcdf files netcdf4 or h5netcdf (default).")
    parser.add_argument('--no_remove_bad', action='store_true',
                        help="Don't remove bad (qual_flags) data.")
    parser.add_argument('-on', '--orbit_number', type=int, nargs='?',
                        required=False, default=99999,
                        help="Orbit number (default is 99999).")
    options = parser.parse_args()
    process_one_file(options.file, options.out_dir, start_line=options.start_line,
                     end_line=options.end_line, engine=options.nc_engine,
                     remove_broken=not options.no_remove_bad,
                     orbit_n=options.orbit_number)
