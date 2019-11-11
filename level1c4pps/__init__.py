#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 level1c4pps developers
#
# This file is part of level1c4pps.
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

#   Adam.Dybbroe <adam.dybbroe@smhi.se>
#   Nina Hakansson <nina.hakansson@smhi.se>

"""Package Initializer for level1c4pps."""

import numpy as np
import xarray as xr
from datetime import datetime
import os
import logging

logger = logging.getLogger('level1c4pps')
from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass


def make_azidiff_angle(sata, suna):
    """Calculate azimuth difference angle."""
    daz = abs(sata - suna)
    daz = daz % 360
    if isinstance(daz, np.ndarray):
        daz[daz > 180] = 360 - daz[daz > 180]
        return daz
    elif isinstance(daz, xr.DataArray):
        return daz.where(daz < 180, 360 - daz)
    else:
        raise ValueError("Azimuth difference is neither a Numpy nor an Xarray object! Type = %s", type(daz))


def dt64_to_datetime(dt64):
    """Conversion of numpy.datetime64 to datetime objects."""
    # https://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64/46921593#46921593
    if type(dt64) == np.datetime64:
        unix_epoch = np.datetime64(0, 's')
        one_second = np.timedelta64(1, 's')
        seconds_since_epoch = (dt64 - unix_epoch) / one_second
        dt = datetime.utcfromtimestamp(seconds_since_epoch)
        return dt
    return dt64

def get_encoding(scene, bandnames, pps_tagnames, angle_names, chunks=None):
    """Get netcdf encoding for all datasets."""
    encoding = {}

    # Bands
    for band in bandnames:
        idtag = pps_tagnames[band]
        try:
            name = scene[band].attrs['name']
        except KeyError:
            logger.debug("No band named %s", band)
            continue
        if 'tb' in idtag:
            encoding[name] = {'dtype': 'int16',
                              'scale_factor': 0.01,
                              '_FillValue': -32767,
                              'zlib': True,
                              'complevel': 4,
                              'add_offset': 273.15}
        else:
            encoding[name] = {'dtype': 'int16',
                              'scale_factor': 0.01,
                              'zlib': True,
                              'complevel': 4,
                              '_FillValue': -32767,
                              'add_offset': 0.0}
        if chunks is not None:
            encoding[name]['chunksizes'] = chunks 

    # Angles and lat/lon
    for name in angle_names:
        encoding[name] = {
            'dtype': 'int16',
            'scale_factor': 0.01,
            'zlib': True,
            'complevel': 4,
            '_FillValue': -32767,
            'add_offset': 0.0}
        if chunks is not None:
            encoding[name]['chunksizes'] = chunks 
    # lat/lon
    for name in ['lon', 'lat']:
        encoding[name] = {'dtype': 'float32',
                          'zlib': True,
                          'complevel': 4,
                          '_FillValue': -999.0}
        if chunks is not None:
            encoding[name]['chunksizes'] = (chunks[1], chunks[2])
    # pygac    
    if 'qual_flags' in scene:
        encoding['qual_flags'] = {'dtype': 'int16', 'zlib': True,
                                  'complevel': 4, '_FillValue': -32001.0}
    if 'scanline_timestamps' in scene:
        encoding['scanline_timestamps'] = {'dtype': 'int64', 'zlib': True,
                                           'complevel': 4, '_FillValue': -1.0}
    return encoding


def compose_filename(scene, out_path, instrument, channel=None):
    """Compose output filename."""
    start_time = scene.attrs['start_time']
    end_time = scene.attrs['end_time']
    if channel is not None:
        start_time = channel.attrs['start_time']
        end_time = channel.attrs['end_time']  
    platform_name = scene.attrs['platform_name']
    orbit_number = int(scene.attrs['orbit_number'])
    filename = os.path.join(
        out_path,
        "S_NWC_{:s}_{:s}_{:05d}_{:s}Z_{:s}Z.nc".format(
            instrument,
            platform_name.lower().replace('-', ''),
            orbit_number,
            datetime.strftime(dt64_to_datetime(start_time), '%Y%m%dT%H%M%S%f')[:-5],

            datetime.strftime(dt64_to_datetime(end_time), '%Y%m%dT%H%M%S%f')[:-5]))
    return filename
