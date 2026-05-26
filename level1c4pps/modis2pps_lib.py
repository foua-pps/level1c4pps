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

"""Functions to convert MODIS level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import time
from satpy.scene import Scene
from level1c4pps import (compose_filename,
                         set_header_and_band_attrs_defaults,
                         rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         convert_angles,
                         save_data,
                         log_time,
                         get_refl_bands,
                         get_band_names,
                         apply_sunz_correction)

import logging
from satpy.utils import debug_on
debug_on()

# Example:
# MYD03.A2014272.0215.006.2014274124038.hdf
# MYD021KM.A2014272.0215.006.2014274142955.hdf


logger = logging.getLogger('modis2pps')

GEOLOCATION_NAMES = [  # additional variables to load
    'satellite_zenith_angle',
    'solar_zenith_angle',
    'satellite_azimuth_angle',
    'solar_azimuth_angle',
    'latitude',
    'longitude']

PPS_TAGNAMES = {'1': 'ch_r06',
                '2': 'ch_r09',
                '26': 'ch_r13',
                '6': 'ch_r16',
                '20': 'ch_tb37',
                '29': 'ch_tb85',
                '31': 'ch_tb11',
                '32': 'ch_tb12',
                # Not used yet:
                '7': 'ch_r21',
                '27': 'ch_tb67',
                '28': 'ch_tb73',
                '33': 'ch_tb133',
                '3': 'ch_rxx',
                '4': 'ch_rxx',
                '5': 'ch_rxx',
                '8': 'ch_rxx',
                '9': 'ch_rxx',
                '10': 'ch_rxx',
                '11': 'ch_rxx',
                '12': 'ch_rxx',
                '13hi': 'ch_rxx',
                '13lo': 'ch_rxx',
                '14hi': 'ch_rxx',
                '14lo': 'ch_rxx',
                '15': 'ch_rxx',
                '16': 'ch_rxx',
                '17': 'ch_rxx',
                '18': 'ch_rxx',
                '19': 'ch_rxx',
                '21': 'ch_tbxx',
                '22': 'ch_tbxx',
                '23': 'ch_tbxx',
                '24': 'ch_tbxx',
                '25': 'ch_tbxx',
                '30': 'ch_tbxx',
                '34': 'ch_tbxx',
                '35': 'ch_tbxx',
                '36': 'ch_tbxx'}

refl_bands = get_refl_bands(PPS_TAGNAMES)
ONE_IR_CHANNEL = '31'


def set_header_and_band_attrs(scene, orbit_n=0):
    """Set and delete some attributes."""
    irch = scene[ONE_IR_CHANNEL]
    nimg = set_header_and_band_attrs_defaults(scene, PPS_TAGNAMES, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "modis2pps.py"
    return nimg


def load_data(scene_files, channel_selection):
    """Load data."""
    scene = Scene(
        reader='modis_l1b',
        filenames=scene_files)
    my_bands = get_band_names(PPS_TAGNAMES, channel_selection)
    bands_to_load = my_bands + GEOLOCATION_NAMES
    scene.load(bands_to_load, resolution=1000)
    return scene


def process_one_scene(scene_files, out_path, engine='h5netcdf', channel_selection="default", orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scene = load_data(scene_files, channel_selection)
    irch = scene[ONE_IR_CHANNEL]
    set_header_and_band_attrs(scene, orbit_n=orbit_n)
    rename_latitude_longitude(scene)
    convert_angles(scene, delete_azimuth=True)
    update_angle_attributes(scene, irch)
    apply_sunz_correction(scene, refl_bands)
    header_attrs = get_header_attrs(scene, band=irch, sensor='modis')
    filename = compose_filename(scene, out_path, instrument='modis', band=irch)
    save_data(scene, filename, header_attrs, engine)
    log_time(filename, tic)
    return filename
