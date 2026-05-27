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

"""Utilities to convert AVHRR GAC formattet data to PPS level-1c format."""


import time
import xarray as xr
import dask.array as da
from datetime import datetime
from satpy.scene import Scene
import pygac  # testing that pygac is available # noqa: F401
from level1c4pps import (log_time,
                         save_data,
                         check_file_exists,
                         compose_filename,
                         set_header_and_band_attrs_defaults,
                         rename_latitude_longitude,
                         update_angle_attributes,
                         get_header_attrs,
                         get_refl_bands,
                         convert_angles)
import logging

logger = logging.getLogger('gac2pps')
GEOLOCATION_NAMES = [  # additional variables to load
    'latitude',
    'longitude',
    'qual_flags',
    'sensor_zenith_angle',
    'solar_zenith_angle',
    'solar_azimuth_angle',
    'sensor_azimuth_angle',
    'sun_sensor_azimuth_difference_angle']
PPS_TAGS = {"1": "ch_r06",
            "2": "ch_r09",
            "3a": "ch_r16",
            "3": "ch_tb37",
            "3b": "ch_tb37",
            "4": "ch_tb11",
            "5": "ch_tb12"}
INSTRUMENTS = {'tirosn': 'avhrr',
               'noaa6': 'avhrr',
               'noaa7': 'avhrr/2',
               'noaa8': 'avhrr',
               'noaa9': 'avhrr/2',
               'noaa10': 'avhrr',
               'noaa11': 'avhrr/2',
               'noaa12': 'avhrr/2',
               'noaa14': 'avhrr/2',
               'noaa15': 'avhrr/3',
               'noaa16': 'avhrr/3',
               'noaa17': 'avhrr/3',
               'noaa18': 'avhrr/3',
               'noaa19': 'avhrr/3'}
refl_bands = get_refl_bands(PPS_TAGS)
band_names = sorted(list(PPS_TAGS.keys()))
ONE_IR_CHANNEL = '4'


def update_ancilliary_datasets(scene):
    """Rename, delete and add some datasets and attributes."""
    irch = scene[ONE_IR_CHANNEL]
    scene['scanline_timestamps'] = xr.DataArray(da.from_array(scene['qual_flags'].coords['acq_time']),
                                                dims=['y'], coords={'y': scene['qual_flags']['y']})
    scene['scanline_timestamps'].attrs['name'] = 'scanline_timestamps'
    # Update qual_flags attrs
    scene['qual_flags'].attrs['id_tag'] = 'qual_flags'
    scene['qual_flags'].attrs['long_name'] = 'pygac quality flags'
    scene['qual_flags'].coords['time'] = irch.attrs['start_time']
    del scene['qual_flags'].coords['acq_time']


def set_header_and_band_attrs(scene, orbit_n=99999):
    """Set and delete some attributes."""
    irch = scene[ONE_IR_CHANNEL]
    nimg = set_header_and_band_attrs_defaults(scene, PPS_TAGS, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "gac2pps.py"
    scene.attrs['is_gac'] = 'True'
    for band in band_names:
        if band not in scene:
            continue
        if band in refl_bands:
            # For GAC data sun_earth_distance_correction is applied always!
            # The sun_earth_distance_correction_factor is not provided by pygac <= 1.2.1 / satpy <= 0.18.1
            scene[band].attrs['sun_earth_distance_correction_applied'] = 'True'
    return nimg


def load_data(gac_file, reader_kwargs):
    """Load data with satpy."""
    if reader_kwargs is None:
        reader_kwargs = {}
    if 'tle_dir' not in reader_kwargs:
        from pygac.configuration import get_config
        conf = get_config()
        tle_dir = conf.get('tle', 'tledir', raw=True)
        tle_name = conf.get('tle', 'tlename', raw=True)
        reader_kwargs['tle_dir'] = tle_dir
        reader_kwargs['tle_name'] = tle_name
    scene = Scene(reader='avhrr_l1b_gaclac',
                  filenames=[gac_file], reader_kwargs=reader_kwargs)
    # Loading all at once sometimes fails with newer satpy, so start with band_names ...
    scene.load(band_names)
    scene.load(GEOLOCATION_NAMES)
    return scene


def process_one_file(gac_file, out_path='.', reader_kwargs=None, engine='h5netcdf', orbit_n=99999):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    check_file_exists(gac_file)
    scene = load_data(gac_file, reader_kwargs)
    irch = scene[ONE_IR_CHANNEL]
    set_header_and_band_attrs(scene, orbit_n=orbit_n)
    rename_latitude_longitude(scene)
    convert_angles(scene)
    update_angle_attributes(scene, irch)
    # Handle gac specific datasets qual_flags and scanline_timestamps
    update_ancilliary_datasets(scene)
    filename = compose_filename(scene, out_path, instrument='avhrr', band=irch)
    header_attrs = get_header_attrs(scene, band=irch, sensor='avhrr')
    save_data(scene, filename, header_attrs, engine)
    log_time(filename, tic)
    return filename
