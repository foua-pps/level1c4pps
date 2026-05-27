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

"""Functions to convert AVHRR AAPP or EPS l1b data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import time
from satpy.scene import Scene
from level1c4pps import (compose_filename,
                         save_data,
                         log_time,
                         rename_latitude_longitude,
                         set_header_and_band_attrs_defaults,
                         update_angle_attributes,
                         get_header_attrs,
                         convert_angles,
                         get_refl_bands,
                         apply_sunz_correction)

import logging

# Example:
# hrpt_noaa19_20141210_1056_30086.l1b

logger = logging.getLogger('avhrr2pps')

GEOLOCATION_NAMES_EPS = [  # additional variables to load
    'satellite_zenith_angle',
    'solar_zenith_angle',
    'satellite_azimuth_angle',
    'solar_azimuth_angle',
    'latitude',
    'longitude']
GEOLOCATION_NAMES_AAPP = [  # additional variables to load
    'sensor_zenith_angle',
    'solar_zenith_angle',
    'sun_sensor_azimuth_difference_angle',
    'latitude',
    'longitude']

PPS_TAGS = {'1': 'ch_r06',
            '2': 'ch_r09',
            '3a': 'ch_r16',
            '3b': 'ch_tb37',
            '4': 'ch_tb11',
            '5': 'ch_tb12'}

refl_bands = get_refl_bands(PPS_TAGS)
band_names = sorted(list(PPS_TAGS.keys()))
ONE_IR_CHANNEL = '4'


def set_header_and_band_attrs(scene, orbit_n=0):
    """Set and delete some attributes."""
    irch = scene[ONE_IR_CHANNEL]
    nimg = set_header_and_band_attrs_defaults(scene, PPS_TAGS, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "avhrr2pps.py"
    return nimg


def check_broken_data(scene):
    """Set bad data to nodata."""
    import numpy as np
    lat = scene['latitude']
    # If we have data in line 2 it is ok
    if (lat[1, :].values > 0).any():
        return
    lon = scene['longitude']
    part_of_lat_that_is_zero = np.sum(lat.values.ravel() == 0) / len(lat.values.ravel())
    part_of_lon_that_is_zero = np.sum(lon.values.ravel() == 0) / len(lon.values.ravel())
    if (part_of_lat_that_is_zero > 0.90 and part_of_lon_that_is_zero > 0.90):
        raise ValueError(
            'More than 90% of data have (lat, lon) equal to (0,0).\n'
            'Most likely this is an older format of hrpt not yet suported by satpy.\n'
            'Stopping. File will not be written.')


def load_data(scene_files):
    """Load data with satpy."""
    if 'AVHR_xxx' in scene_files[0]:
        avhrr_reader = 'avhrr_l1b_eps'
        angles = GEOLOCATION_NAMES_EPS
    else:
        avhrr_reader = 'avhrr_l1b_aapp'
        angles = GEOLOCATION_NAMES_AAPP
    scene = Scene(
        reader=avhrr_reader,
        filenames=scene_files)
    bands_to_load = band_names + angles
    scene.load(bands_to_load)
    return scene


def process_one_scene(scene_files, out_path, engine='h5netcdf', orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scene = load_data(scene_files)
    irch = scene[ONE_IR_CHANNEL]
    # Check if we have old hrpt format with data only every 20th line
    check_broken_data(scene)
    set_header_and_band_attrs(scene, orbit_n=orbit_n)
    rename_latitude_longitude(scene)
    convert_angles(scene, delete_azimuth=True)
    update_angle_attributes(scene, irch)
    apply_sunz_correction(scene, refl_bands)
    header_attrs = get_header_attrs(scene, band=irch, sensor='avhrr')
    filename = compose_filename(scene, out_path, instrument='avhrr', band=irch)
    save_data(scene, filename, header_attrs, engine)
    log_time(filename, tic)
    return filename
