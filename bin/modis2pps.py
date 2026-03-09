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

"""Script to convert MODIS level-1 to PPS level-1c format using Pytroll/Satpy."""

import argparse
from level1c4pps.modis2pps_lib import process_one_scene


if __name__ == "__main__":
    """ Create PPS-format level1c data
    From a list of MERSI-2 level-1 files create a NWCSAF/PPS formatet level1c file for pps.
    """
    parser = argparse.ArgumentParser(
        description=('Script to produce a PPS-level1c file for a MODIS level-1 scene'))
    parser.add_argument('files', metavar='fileN', type=str, nargs='+',
                        help='List of MODIS files to process')
    parser.add_argument('-o', '--out_dir', type=str, nargs='?',
                        required=False, default='.',
                        help="Output directory where to store the level1c file")
    parser.add_argument('-ne', '--nc_engine', type=str, nargs='?',
                        required=False, default='h5netcdf',
                        help="Engine for saving netcdf files netcdf4 or h5netcdf (default).")
    parser.add_argument('-all_ch', '--all_channels', action='store_true',
                        help="Save all 36 channels to level1c4pps file.")
    parser.add_argument('-pps_ch', '--pps_channels', action='store_true',
                        help="Save only the necessary (for PPS) channels to level1c4pps file.")
    parser.add_argument('-on', '--orbit_number', type=int, nargs='?',
                        required=False, default=0,
                        help="Orbit number (default is 00000).")

    options = parser.parse_args()
    process_one_scene(options.files, options.out_dir, engine=options.nc_engine,
                      all_channels=options.all_channels, pps_channels=options.pps_channels,
                      orbit_n=options.orbit_number)
