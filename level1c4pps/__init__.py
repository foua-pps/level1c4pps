#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 level1c4pps developers
#
# This file is part of level1c4pps.
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

#   Adam.Dybbroe <adam.dybbroe@smhi.se>
#   Nina Hakansson <nina.hakansson@smhi.se>

"""Package Initializer for level1c4pps."""

from pkg_resources import get_distribution, DistributionNotFound
import numpy as np
import xarray as xr
from datetime import datetime
import os
import logging

logger = logging.getLogger('level1c4pps')
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

PPS_ANGLE_TAGS = ['sunzenith', 'satzenith', 'azimuthdiff', 'sunazimuth', 'satazimuth']
ANGLE_ATTRIBUTES = {
    'long_name': {
        'sunzenith': 'sun zenith angle',
        'satzenith': 'satellite zenith angle',
        'azimuthdiff': 'absolute azimuth difference angle',
        'sunazimuth': 'sun azimuth angle degree clockwise from north',
        'satazimuth': 'satellite azimuth angle degree clockwise from north',
    },
    'valid_range': {
        'sunzenith': [0, 18000],
        'satzenith': [0, 9000],
        'azimuthdiff': [0, 18000],
        'sunazimuth': [-18000, 18000],
        'satazimuth': [-18000, 18000],
    },
    'mersi2_file_key':  {
        'sunzenith': 'Geolocation/SolarZenithAngle',
        'satzenith': 'Geolocation/SensorZenithAngle',
        'azimuthdiff': 'Geolocation/SensorSolarAzimuthDifference',
    },
    'standard_name':  {
        'sunzenith': 'solar_zenith_angle',
        'satzenith': 'sensor_zenith_angle',  # platform in ppsv2018
        'azimuthdiff': 'absolute_angle_of_rotation_from_solar_azimuth_to_platform_azimuth',
        'sunazimuth': 'solar_azimuth_angle',
        'satazimuth': 'sensor_azimuth_angle',  # plaform in ppsv2018
    }
}


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


def get_encoding(scene, bandnames, pps_tagnames, chunks=None):
    """Get netcdf encoding for all datasets."""
    encoding = {}
    for dataset in scene.keys():
        try:
            name, enc = get_band_encoding(scene[dataset.name], bandnames, pps_tagnames,
                                          chunks=chunks)
        except ValueError:
            continue
        encoding[name] = enc
    return encoding


def get_band_encoding(dataset, bandnames, pps_tagnames, chunks=None):
    """Get netcdf encoding for a datasets."""
    name = dataset.attrs['name']
    id_tag = dataset.attrs.get('id_tag', None)
    if id_tag is not None:
        if id_tag.startswith('ch_tb'):
            # IR channel
            enc = {'dtype': 'int16',
                   'scale_factor': 0.01,
                   '_FillValue': -32767,
                   'zlib': True,
                   'complevel': 4,
                   'add_offset': 273.15}
        if id_tag.startswith('ch_r'):
            # Refl channel
            enc = {'dtype': 'int16',
                   'scale_factor': 0.01,
                   'zlib': True,
                   'complevel': 4,
                   '_FillValue': -32767,
                   'add_offset': 0.0}
        if id_tag in PPS_ANGLE_TAGS:
            # Angle
            enc = {
                'dtype': 'int16',
                'scale_factor': 0.01,
                'zlib': True,
                'complevel': 4,
                '_FillValue': -32767,
                'add_offset': 0.0}
        if chunks is not None:
            enc['chunksizes'] = chunks
    elif name in ['lon', 'lat']:
        # Lat/Lon
        enc = {'dtype': 'float32',
               'zlib': True,
               'complevel': 4,
               '_FillValue': -999.0}
        if chunks is not None:
            enc['chunksizes'] = (chunks[1], chunks[2])
    elif name in ['qual_flags']:
        # pygac qual flags
        enc = {'dtype': 'int16', 'zlib': True,
               'complevel': 4, '_FillValue': -32001.0}
    elif name in ['scanline_timestamps']:
        # pygac scanline_timestamps
        enc = {'dtype': 'int64', 'zlib': True,
               'complevel': 4, '_FillValue': -1.0}
    else:
        raise ValueError('Unsupported band: {}'.format(name))
    return name, enc


def rename_latitude_longitude(scene):
    """Rename latitude longitude to lat lon."""
    scene['lat'] = scene['latitude']
    scene['lon'] = scene['longitude']
    del scene['latitude']
    del scene['longitude']
    # Update attributes
    scene['lat'].attrs['long_name'] = 'latitude coordinate'
    scene['lon'].attrs['long_name'] = 'longitude coordinate'
    scene['lat'].attrs['name'] = 'lat'
    scene['lon'].attrs['name'] = 'lon'
    try:
        del scene['lat'].coords['acq_time']
        del scene['lon'].coords['acq_time']
    except KeyError:
        pass


def update_angle_attributes(scene, band):
    """Set and delete angle attributes."""
    for angle in PPS_ANGLE_TAGS:
        if angle not in scene.keys() and angle in ['sunazimuth', 'satazimuth']:
            # azimuth angles not always there
            continue
        scene[angle].attrs['id_tag'] = angle
        scene[angle].attrs['name'] = angle
        scene[angle].attrs['coordinates'] = 'lon lat'
        scene[angle].attrs['units'] = 'degree'
        scene[angle].attrs['long_name'] = ANGLE_ATTRIBUTES['long_name'][angle]
        scene[angle].attrs['valid_range'] = ANGLE_ATTRIBUTES['valid_range'][angle]
        scene[angle].attrs['standard_name'] = ANGLE_ATTRIBUTES['standard_name'][angle]
        scene[angle].coords['time'] = band.attrs["start_time"]
        for attr in ["start_time", "end_time"]:
            scene[angle].attrs[attr] = band.attrs[attr]
        # delete some attributes
        try:
            del scene[angle].attrs['area']
        except KeyError:
            pass
        # delete some coords
        try:
            del scene[angle].coords['acq_time']
        except KeyError:
            pass


def compose_filename(scene, out_path, instrument, band=None):
    """Compose output filename.

    As default use the start and end time of the scene.
    For SEVIRI this is the nominal timestamp of the scan (as in the HRIT files).
    If a scene band is supplied use that for start/end time.

    Args:
        scene: satpy scene
        outpath: output directory (string)
        instrument: lower case instrument (string)
        band: use start and end time from band if supplied

    """
    start_time = scene.attrs['start_time']
    end_time = scene.attrs['end_time']
    if band is not None:
        start_time = band.attrs['start_time']
        end_time = band.attrs['end_time']
    platform_name = scene.attrs['platform']
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


def get_header_attrs(scene, band, sensor='avhrr'):
    """Get global netcdf attributes."""
    header_attrs = scene.attrs.copy()
    header_attrs['start_time'] = datetime.strftime(dt64_to_datetime(band.attrs['start_time']),
                                                   "%Y-%m-%d %H:%M:%S")
    header_attrs['end_time'] = datetime.strftime(dt64_to_datetime(band.attrs['end_time']),
                                                 "%Y-%m-%d %H:%M:%S")
    header_attrs['sensor'] = sensor
    return header_attrs
