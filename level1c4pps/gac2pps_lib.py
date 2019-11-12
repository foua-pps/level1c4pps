#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 level1c4pps developers
#
# This file is part of level1c4pps
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
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
from level1c4pps import dt64_to_datetime, get_encoding, compose_filename
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


def get_encoding_gac(scene, angle_names):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        angle_names,
                        chunks=None)


def update_ancilliary_datasets(scene, image_num):
    """Rename, delete and add some datasets and attributes."""
    irch = scene['4']

    # Rename latitude longitude
    scene['lat'] = scene['latitude']
    scene['lon'] = scene['longitude']
    del scene['latitude']
    del scene['longitude']
    # Update attributes
    scene['lat'].attrs['long_name'] = 'latitude coordinate'
    scene['lon'].attrs['long_name'] = 'longitude coordinate'
    del scene['lat'].coords['acq_time']
    del scene['lon'].coords['acq_time']

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

    # Rename angle datasets
    scene['sunzenith'] = scene['solar_zenith_angle']
    scene['satzenith'] = scene['sensor_zenith_angle']
    # PPS need absolute azimuth angle differences
    scene['azimuthdiff'] = abs(scene['sun_sensor_azimuth_difference_angle'])
    scene['azimuthdiff'].attrs = scene['sun_sensor_azimuth_difference_angle'].attrs
    scene['sunazimuth'] = scene['solar_azimuth_angle']
    scene['satazimuth'] = scene['sensor_azimuth_angle']
    del scene['solar_zenith_angle']
    del scene['sensor_zenith_angle']
    del scene['sun_sensor_azimuth_difference_angle']
    del scene['solar_azimuth_angle']
    del scene['sensor_azimuth_angle']

    # Set angle specific attributes
    scene['sunzenith'].attrs['long_name'] = 'sun zenith angle'
    scene['sunzenith'].attrs['valid_range'] = [0, 18000]
    scene['satzenith'].attrs['long_name'] = 'satellite zenith angle'
    scene['satzenith'].attrs['valid_range'] = [0, 9000]
    scene['azimuthdiff'].attrs['long_name'] = 'absolute azimuth difference angle'
    scene['azimuthdiff'].attrs['valid_range'] = [0, 18000]
    scene['sunazimuth'].attrs['long_name'] = 'sun azimuth angle degree clockwise from north'
    scene['sunazimuth'].attrs['valid_range'] = [-18000, 18000]
    scene['satazimuth'].attrs['long_name'] = 'satellite azimuth angle degree clockwise from north'
    scene['satazimuth'].attrs['valid_range'] = [-18000, 18000]

    angle_names = []
    for image_num, angle in enumerate(['sunzenith', 'satzenith', 'azimuthdiff', 'sunazimuth', 'satazimuth']):
        scene[angle].attrs['id_tag'] = angle
        scene[angle].attrs['name'] = "image{:d}".format(image_num)
        angle_names.append("image{:d}".format(image_num))
        image_num += 1
        scene[angle].attrs['coordinates'] = 'lon lat'
        scene[angle].coords['time'] = irch.attrs['start_time']
        # delete some attributes
        del scene[angle].coords['acq_time']
        del scene[angle].attrs['area']
    return angle_names


def set_header_and_band_attrs(scene):
    """Set and delete some attributes."""

    # Set some header attributes:
    scene.attrs['instrument'] = "AVHRR"
    scene.attrs['source'] = "gac2pps.py"
    nowutc = datetime.utcnow()
    scene.attrs['date_created'] = nowutc.strftime("%Y-%m-%dT%H:%M:%SZ")
    irch = scene['4']
    scene.attrs['platform'] = irch.attrs['platform_name']
    scene.attrs['platform_name'] = irch.attrs['platform_name']
    scene.attrs['orbit_number'] = '{:05d}'.format(irch.attrs['orbit_number'])
    scene.attrs['orbit'] = scene.attrs['orbit_number']

    # bands
    image_num = 0 # name of first dataset is image0
    for band in enumerate(BANDNAMES):
        try:
            idtag = PPS_TAGNAMES.get(band, band)
            scene[band].attrs['id_tag'] = idtag
            scene[band].attrs['description'] = 'AVHRR ' + str(band)
            scene[band].attrs['sun_earth_distance_correction_applied'] = 'False'
            scene[band].attrs['sun_earth_distance_correction_factor'] = 1.0
            scene[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
            scene[band].attrs['name'] = "image{:d}".format(image_num)
            image_num += 1
            scene[band].attrs['coordinates'] = 'lon lat'
            del scene[band].attrs['area']
        except KeyError:
            continue
    return image_num


def get_header_attrs(scene):
    """Get global netcdf attributes."""
    header_attrs = scene.attrs.copy()
    irch = scene['4']
    header_attrs['start_time'] = datetime.strftime(dt64_to_datetime(irch.attrs['start_time']),
                                                   "%Y-%m-%d %H:%M:%S")
    header_attrs['end_time'] = datetime.strftime(dt64_to_datetime(irch.attrs['end_time']),
                                                 "%Y-%m-%d %H:%M:%S")
    header_attrs['sensor'] = 'avhrr'
    return header_attrs


def process_one_file(gac_file, out_path='.', reader_kwargs=None):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(reader='avhrr_l1b_gaclac',
                 filenames=[gac_file], reader_kwargs=reader_kwargs)
    if 'avhrr-3' in scn_.attrs['sensor']:
        sensor = 'avhrr'
        scn_.load(BANDNAMES + ['latitude', 'longitude',
                               'qual_flags',
                               'sensor_zenith_angle', 'solar_zenith_angle',
                               'solar_azimuth_angle', 'sensor_azimuth_angle',
                               'sun_sensor_azimuth_difference_angle'])

    # Set header and band attributes
    image_num = set_header_and_band_attrs(scn_)

    # Rename and set att header and band attributes
    angle_names = update_ancilliary_datasets(scn_, image_num)

    filename = compose_filename(scn_, out_path, instrument='avhrr', band=scn_['4'])
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_),
                       engine='netcdf4',
                       flatten_attrs=True,
                       encoding=get_encoding_gac(scn_, angle_names))

    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
