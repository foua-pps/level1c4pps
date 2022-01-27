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

"""Script to convert METIMAGE level-1 to PPS level-1c format using Pytroll/Satpy."""

import argparse
from level1c4pps.metimage2pps_lib import process_one_scene


if __name__ == "__main__":
    """ Create PPS-format level1c data
    From a list of MERSI-2 level-1 files create a NWCSAF/PPS formatet level1c file for pps.
    """
    parser = argparse.ArgumentParser(
        description=('Script to produce a PPS-level1c file for a METIMAGE level-1 scene'))
    parser.add_argument('files', metavar='fileN', type=str, nargs='+',
                        help='List of METIMAGE files to process')
    parser.add_argument('-o', '--out_dir', type=str, nargs='?',
                        required=False, default='.',
                        help="Output directory where to store the level1c file")
    parser.add_argument('-all_ch', '--all_channels', action='store_true',
                        help="Save all 20 channels to level1c4pps file.")
    parser.add_argument('-pps_ch', '--pps_channels', action='store_true',
                        help="Save only the necessary (for PPS) channels to level1c4pps file.")
    parser.add_argument('-ne', '--nc_engine', type=str, nargs='?',
                        required=False, default='h5netcdf',
                        help="Engine for saving netcdf files netcdf4 or h5netcdf (default).")
    parser.add_argument('-on', '--orbit_number', type=int, nargs='?',
                        required=False, default=0,
                        help="Orbit number (default is 00000).")
    options = parser.parse_args()
    process_one_scene(options.files, options.out_dir,
                      engine=options.nc_engine,
                      all_channels=options.all_channels,
                      pps_channels=options.pps_channels,
                      orbit_n=options.orbit_number)
