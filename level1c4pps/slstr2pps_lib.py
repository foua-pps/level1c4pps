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

"""Functions to convert MERSI-2 level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import time
import satpy
from satpy.scene import Scene
from level1c4pps import (save_data,
                         log_time,
                         compose_filename,
                         set_header_and_band_attrs_defaults,
                         rename_latitude_longitude,
                         update_angle_attributes,
                         get_band_names,
                         get_refl_bands,
                         get_header_attrs,
                         convert_angles)
import pyspectral  # testing that pyspectral is available # noqa: F401
import logging
from packaging.version import Version

# Example:
# S3A_SL_1_RBT____20200911T163815_20200911T164115_20200912T215136_0179_062_340_0360_LN2_O_NT_004.SEN3#

logger = logging.getLogger('slstr2pps')

# slstr2pps.py -o temp/ S3A_SL_1_RBT____20200911T163815_*.SEN3/*n.nc
# ['F1_in', 'F2_in', 'S1_an', 'S2_an', 'S3_an', 'S4_an', 'S4_bn', 'S5_an', 'S5_bn',
# 'S6_an', 'S6_bn', 'S7_in', 'S8_in', 'S9_in',
# 'bayes_an', 'bayes_bn', 'bayes_in', 'cloud_an', 'cloud_bn', 'cloud_in',
# 'confidence_an', 'confidence_bn', 'confidence_in',
# 'latitude_an', 'latitude_ao', 'latitude_bn', 'latitude_bo', 'latitude_in', 'latitude_io',
# 'longitude_an', 'longitude_ao', 'longitude_bn', 'longitude_bo', 'longitude_in', 'longitude_io',
# 'pointing_an', 'pointing_bn', 'pointing_in',
# 'satellite_azimuth_angle_n', 'satellite_azimuth_angle_o', 'satellite_zenith_angle_n', 'satellite_zenith_angle_o',
# 'solar_azimuth_angle_n', 'solar_azimuth_angle_o', 'solar_zenith_angle_n', 'solar_zenith_angle_o']


PPS_TAGNAMES = {'S2': 'ch_r06',  # or S1
                'S3': 'ch_r09',
                'S4': 'ch_r13',
                'S5': 'ch_r16',
                'S6': 'ch_r22',
                'S7': 'ch_tb37',
                'S8': 'ch_tb11',
                'S9': 'ch_tb12',
                # Not yet in pps:
                'S1': 'ch_rxx',
                'F1': 'ch_tbxx',
                'F2': 'ch_tbxx'}
refl_bands = get_refl_bands(PPS_TAGNAMES)

GEOLOCATION_NAMES = [  # additional variables to load
    'latitude',
    'longitude',
    'satellite_azimuth_angle',
    'satellite_zenith_angle',
    'solar_azimuth_angle',
    'solar_zenith_angle',
]

ONE_IR_CHANNEL = 'S8'


def set_header_and_band_attrs(scene, orbit_n=0):
    """Set and delete some attributes."""
    irch = scene[ONE_IR_CHANNEL]
    nimg = set_header_and_band_attrs_defaults(scene, PPS_TAGNAMES, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "slstr2pps.py"
    return nimg


def load_data(scene_files, all_channels=False, pps_channels=False):
    """Load data with satpy."""
    scene = Scene(
        reader='slstr_l1b',
        filenames=scene_files)
    my_bands = get_band_names(PPS_TAGNAMES, all_channels, pps_channels)
    scene.load(my_bands + GEOLOCATION_NAMES)
    # Everything should be on the same grid, to be saved as ppsleve1c
    scene = scene.resample(resampler="native")


def process_one_scene(scene_files, out_path, engine='h5netcdf',
                      all_channels=False, pps_channels=False, orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scene = load_data(scene_files, all_channels=all_channels, pps_channels=pps_channels)
    irch = scene[ONE_IR_CHANNEL]
    set_header_and_band_attrs(scene, orbit_n=orbit_n)
    rename_latitude_longitude(scene)
    convert_angles(scene, delete_azimuth=True)
    update_angle_attributes(scene, irch)
    header_attrs = get_header_attrs(scene, band=irch, sensor='slstr')
    filename = compose_filename(scene, out_path, instrument='slstr', band=irch)
    save_data(scene, filename, header_attrs, engine)
    log_time(filename, tic)
    return filename
