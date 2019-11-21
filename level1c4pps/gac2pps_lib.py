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

"""Utilities to convert AVHRR GAC formattet data to PPS level-1c format."""


import os
import time
import xarray as xr
import dask.array as da
import numpy as np
from datetime import datetime
from satpy.scene import Scene
import pygac  # testing that pygac is available # noqa: F401
from level1c4pps import (get_encoding, compose_filename,
                         rename_latitude_longitude, update_angle_attributes,
                         get_header_attrs)
import logging

logger = logging.getLogger('gac2pps')

BANDNAMES = ['1', '2', '3', '3a', '3b', '4', '5']

PPS_TAGNAMES = {"1": "ch_r06",
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
SATPY_ANGLE_NAMES = {
    'sunzenith': 'solar_zenith_angle',
    'satzenith': 'sensor_zenith_angle',
    'azimuthdiff': 'sun_sensor_azimuth_difference_angle',
    'sunazimuth': 'solar_azimuth_angle',
    'satazimuth': 'sensor_azimuth_angle'}


def get_encoding_gac(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def update_ancilliary_datasets(scene):
    """Rename, delete and add some datasets and attributes."""
    irch = scene['4']

    # Create new data set scanline timestamps
    first_jan_1970 = np.array([datetime(1970, 1, 1, 0, 0, 0)]).astype('datetime64[ns]')
    scanline_timestamps = np.array(scene['qual_flags'].coords['acq_time'] -
                                   first_jan_1970).astype(dtype='timedelta64[ms]').astype(np.float64)
    scene['scanline_timestamps'] = xr.DataArray(da.from_array(scanline_timestamps),
                                                dims=['y'], coords={'y': scene['qual_flags']['y']})
    scene['scanline_timestamps'].attrs['units'] = 'Milliseconds since 1970-01-01 00:00:00 UTC'

    # Update qual_flags attrs
    scene['qual_flags'].attrs['id_tag'] = 'qual_flags'
    scene['qual_flags'].attrs['long_name'] = 'pygac quality flags'
    scene['qual_flags'].coords['time'] = irch.attrs['start_time']
    del scene['qual_flags'].coords['acq_time']


def convert_angles(scene, satpy_angle_names):
    """Convert angles to pps format."""
    for angle in ['sunzenith', 'satzenith', 'sunazimuth', 'satazimuth']:
        scene[angle] = scene[satpy_angle_names[angle]]
        del scene[satpy_angle_names[angle]]
    angle = 'azimuthdiff'
    scene[angle] = abs(scene[satpy_angle_names[angle]])
    scene[angle].attrs = scene[satpy_angle_names[angle]].attrs
    del scene[satpy_angle_names[angle]]


def set_header_and_band_attrs(scene):
    """Set and delete some attributes."""
    scene.attrs['instrument'] = "AVHRR"
    scene.attrs['source'] = "gac2pps.py"
    scene.attrs['history'] = "Created by level1c4pps."  # history attr missing in satpy=0.18.2
    nowutc = datetime.utcnow()
    scene.attrs['date_created'] = nowutc.strftime("%Y-%m-%dT%H:%M:%SZ")
    irch = scene['4']
    scene.attrs['platform'] = irch.attrs['platform_name']
    scene.attrs['platform_name'] = irch.attrs['platform_name']
    scene.attrs['orbit_number'] = irch.attrs['orbit_number']
    scene.attrs['orbit'] = scene.attrs['orbit_number']

    # bands
    image_num = 0  # name of first dataset is image0
    for band in BANDNAMES:
        try:
            idtag = PPS_TAGNAMES.get(band, band)
            scene[band].attrs['id_tag'] = idtag
            scene[band].attrs['description'] = 'AVHRR ' + str(band)
            scene[band].attrs['sun_earth_distance_correction_applied'] = 'True'
            scene[band].attrs['sun_earth_distance_correction_factor'] = 1.0
            scene[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
            scene[band].attrs['name'] = "image{:d}".format(image_num)
            image_num += 1
            scene[band].attrs['coordinates'] = 'lon lat'
            # Add time coordinate. To make cfwriter aware that we want 3D data.
            scene[band].coords['time'] = irch.attrs['start_time']
            del scene[band].coords['acq_time']
            del scene[band].attrs['area']
        except KeyError:
            continue


def process_one_file(gac_file, out_path='.', reader_kwargs=None):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(reader='avhrr_l1b_gaclac',
                 filenames=[gac_file], reader_kwargs=reader_kwargs)

    scn_.load(BANDNAMES + ['latitude', 'longitude',
                           'qual_flags',
                           'sensor_zenith_angle', 'solar_zenith_angle',
                           'solar_azimuth_angle', 'sensor_azimuth_angle',
                           'sun_sensor_azimuth_difference_angle'])
    # one ir channel
    irch = scn_['4']

    # Set header and band attributes
    set_header_and_band_attrs(scn_)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, SATPY_ANGLE_NAMES)
    update_angle_attributes(scn_, irch)

    # Handle gac specific datasets qual_flags and scanline_timestamps
    update_ancilliary_datasets(scn_)

    filename = compose_filename(scn_, out_path, instrument='avhrr', band=irch)
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='avhrr'),
                       engine='netcdf4',
                       flatten_attrs=True,
                       encoding=get_encoding_gac(scn_))

    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
