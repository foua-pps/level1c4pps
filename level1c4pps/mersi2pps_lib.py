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

"""Functions to convert MERSI-2/3 level-1 data to a NWCSAF/PPS level-1c formated netCDF/CF file."""

import logging
import os
import time

import numpy as np
from satpy.scene import Scene

from level1c4pps import (ANGLE_ATTRIBUTES, check_file_exists, compose_filename,
                         convert_angles, get_band_names, get_header_attrs,
                         get_refl_bands, log_time, rename_latitude_longitude,
                         save_data, set_header_and_band_attrs_defaults,
                         update_angle_attributes)

# Example filenames:
# tf2019234102243.FY3D-X_MERSI_GEOQK_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_GEO1K_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_1000M_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_0250M_L1B.HDF
#

logger = logging.getLogger('mersi2pps')

SENSOR = {  # sensor associated to platform
    'FY3D': 'mersi-2',
    'FY3F': 'mersi-3',
    'FY3H': 'mersi-3',
}
SATPY_READER = {  # satpy reader associated to sensor
    'mersi-2': 'mersi2_l1b',
    'mersi-3': 'mersi3_l1b',
}
PPS_TAGS = {  # PPS band name associated to satpy band name
    '3': 'ch_r06',
    '4': 'ch_r09',
    '5': 'ch_r13',
    '6': 'ch_r16',
    '20': 'ch_tb37',
    '22': 'ch_tb73',
    '23': 'ch_tb85',
    '24': 'ch_tb11',
    '25': 'ch_tb12',
}
GEOLOCATION_NAMES = [  # additional variables to load
    'latitude',
    'longitude',
    'satellite_azimuth_angle',
    'satellite_zenith_angle',
    'solar_azimuth_angle',
    'solar_zenith_angle',
]

refl_bands = get_refl_bands(PPS_TAGS)

ONE_IR_CHANNEL = '24'
RESOLUTION = 1000  # [m]
LOW_TB = 1  # [K] very cold brightness temperature


def set_header_and_band_attrs(scene, ir_channel_obj, orbit_n):
    """Set and delete some attributes."""
    set_header_and_band_attrs_defaults(
        scene, PPS_TAGS, ir_channel_obj, orbit_n=orbit_n,
    )
    scene.attrs['source'] = "mersi2pps.py"


def remove_broken_data(scene):
    """Set bad data to nodata."""
    for band in PPS_TAGS:
        if band not in refl_bands and band in scene:
            scene[band].data = np.where(scene[band].data < LOW_TB, np.nan, scene[band].data)


def get_sensor(scene_file):
    """Get sensor associated to the scene file."""
    for platform, sensor in SENSOR.items():
        if platform in scene_file:
            return sensor
    logger.info(f"Failed to determine sensor associated to scene file: {scene_file}")
    return None


def load_data(scene_files, sensor, channel_selection):
    """Load data."""
    reader = SATPY_READER[sensor]
    scene = Scene(reader=reader, filenames=scene_files)
    band_names = get_band_names(PPS_TAGS, channel_selection)
    bands_to_load = band_names + GEOLOCATION_NAMES
    scene.load(bands_to_load, resolution=RESOLUTION)
    return scene


def process_one_scene(scene_files, out_path, channel_selection="default", engine='h5netcdf', orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    check_file_exists(scene_files)
    sensor = get_sensor(os.path.basename(scene_files[0]))
    scene = load_data(scene_files, sensor, channel_selection)
    remove_broken_data(scene)
    ir_channel_obj = scene[ONE_IR_CHANNEL]
    set_header_and_band_attrs(scene, ir_channel_obj, orbit_n)
    rename_latitude_longitude(scene)
    convert_angles(scene, delete_azimuth=True)
    update_angle_attributes(scene, ir_channel_obj)
    for angle in ['sunzenith', 'satzenith', 'azimuthdiff']:
        scene[angle].attrs['file_key'] = ANGLE_ATTRIBUTES['mersi_file_key'][angle]
    header_attrs = get_header_attrs(scene, band=ir_channel_obj, sensor=sensor)
    filename = compose_filename(scene, out_path, instrument=sensor.replace('-', ''), band=band)
    save_data(scene, filename, header_attrs, engine)
    log_time(filename, tic)
    return filename
