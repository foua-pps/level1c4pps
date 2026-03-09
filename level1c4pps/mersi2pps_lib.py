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

"""Functions to convert MERSI-2/3 level-1 data to a NWCSAF/PPS level-1c formated netCDF/CF file."""

import logging
import os
import time

import numpy as np
from satpy.scene import Scene
from level1c4pps import (
    ANGLE_ATTRIBUTES,
    compose_filename,
    convert_angles,
    get_encoding,
    get_header_attrs,
    rename_latitude_longitude,
    set_header_and_band_attrs_defaults,
    update_angle_attributes,
)

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
PPS_BAND_NAME = {  # PPS band name associated to satpy band name
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
REFL_BANDS = ['3', '4', '5', '6']
ONE_IR_CHANNEL = '24'
RESOLUTION = 1000  # [m]
LOW_TB = 1  # [K] very cold brightness temperature


def set_header_and_band_attrs(scene, band, orbit_n):
    """Set and delete some attributes."""
    set_header_and_band_attrs_defaults(
        scene, list(PPS_BAND_NAME), PPS_BAND_NAME, REFL_BANDS, band, orbit_n=orbit_n,
    )
    scene.attrs['source'] = "mersi2pps.py"


def remove_broken_data(scene):
    """Set bad data to nodata."""
    for band in PPS_BAND_NAME:
        if band not in REFL_BANDS and band in scene:
            scene[band].data = np.where(scene[band].data < LOW_TB, np.nan, scene[band].data)


def get_sensor(scene_file):
    """Get sensor associated to the scene file."""
    for platform, sensor in SENSOR.items():
        if platform in scene_file:
            return sensor
    logger.info("Failed to determine sensor associated to scene file: '%s'", scene_file)
    return None


def process_one_scene(scene_files, out_path, engine='h5netcdf', orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    sensor = get_sensor(os.path.basename(scene_files[0]))
    reader = SATPY_READER[sensor]
    scene = Scene(reader=reader, filenames=scene_files)
    band_names = list(PPS_BAND_NAME)
    scene.load(band_names + GEOLOCATION_NAMES, resolution=RESOLUTION)
    remove_broken_data(scene)
    band = scene[ONE_IR_CHANNEL]
    set_header_and_band_attrs(scene, band, orbit_n)
    rename_latitude_longitude(scene)
    convert_angles(scene, delete_azimuth=True)
    update_angle_attributes(scene, band)
    for angle in ['sunzenith', 'satzenith', 'azimuthdiff']:
        scene[angle].attrs['file_key'] = ANGLE_ATTRIBUTES['mersi_file_key'][angle]
    filename = compose_filename(scene, out_path, instrument=sensor.replace('-', ''), band=band)
    scene.save_datasets(
        writer='cf',
        filename=filename,
        header_attrs=get_header_attrs(scene, band=band, sensor=sensor),
        engine=engine,
        include_lonlats=False,
        flatten_attrs=True,
        encoding=get_encoding(scene, band_names, PPS_BAND_NAME, chunks=None),
    )
    print(f"Saved file {os.path.basename(filename)} after {time.time() - tic:3.1f} seconds")
    return filename
