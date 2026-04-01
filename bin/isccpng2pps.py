#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2025 level1c4pps developers
#
# This file is part of level1c4pps
#
# level1c4pps is free software: you can redistribute it and/or modify it
# under the terms of the Apache-2.0 license.
#
#
# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Nina Hakansson <nina.hakansson@smhi.se>
#   Adam.Dybbroe <adam.dybbroe@smhi.se>


"""Script to make ISCCP-NG level1c in PPS-format with pytroll."""

import argparse
from level1c4pps.isccpng2pps_lib import process_one_scene


# -----------------------------------------------------------------------------
# Main:
if __name__ == "__main__":
    """ Create PPS-format level1c data
    From a list of hirt files hrit create a level1c file for pps.
    """
    parser = argparse.ArgumentParser(
        description=('Script to produce a PPS-level1c file for a '
                     'ISCCP-NG scene.'))
    parser.add_argument('files', metavar='fileN', type=str, nargs='+',
                        help='List of ISCCP-NG files to process')
    parser.add_argument('-o', '--out_dir', type=str, nargs='?',
                        required=False, default='.',
                        help="Output directory where to store level1c file, can include time directives (%Y%m%d).")
    parser.add_argument('-ne', '--nc_engine', type=str, nargs='?',
                        required=False, default='h5netcdf',
                        help="Engine for saving netcdf files netcdf4 or h5netcdf (default).")
    parser.add_argument('-on', '--orbit_number', type=int, nargs='?',
                        required=False, default=99999,
                        help="Orbit number (default is 99999).")
    options = parser.parse_args()
    process_one_scene(options.files, options.out_dir,
                      engine=options.nc_engine,
                      orbit_n=options.orbit_number)
