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
import satpy
import level1c4pps

logger = logging.getLogger('level1c4pps')
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

SATPY_ANGLE_NAMES = {
    'solar_zenith': 'sunzenith',  # no _angle
    'solar_zenith_angle': 'sunzenith',
    'solar_azimuth': 'sunazimuth',  # no _angle
    'solar_azimuth_angle': 'sunazimuth',
    'satellite_zenith_angle': 'satzenith',
    'sensor_zenith_angle': 'satzenith',
    'observation_zenith': 'satzenith',
    'satellite_azimuth_angle': 'satazimuth',
    'sensor_azimuth_angle': 'satazimuth',
    'observation_azimuth': 'satazimuth',
    'sun_sensor_azimuth_difference_angle': 'azimuthdiff',
}


def convert_angles(scene, delete_azimuth=False):
    """Convert angles to pps format."""
    for satpy_name in SATPY_ANGLE_NAMES:
        if satpy_name in scene:
            scene[SATPY_ANGLE_NAMES[satpy_name]] = scene[satpy_name]  # Rename angle
            del scene[satpy_name]

    angle = 'azimuthdiff'
    if angle not in scene:
        # Create azimuth diff angle
        scene[angle] = make_azidiff_angle(scene['satazimuth'], scene['sunazimuth'])
        scene[angle].attrs = scene['sunazimuth'].attrs  # Copy sunazimuth attrs
    else:
        # Just apply abs
        scene[angle] = abs(scene[angle])
        scene[angle].attrs = scene[angle].attrs

    if delete_azimuth:
        # PPS does not need azimuth angles
        try:
            del scene['satazimuth']
            del scene['sunazimuth']
        except KeyError:
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
        'sunzenith': np.array([0, 18000], dtype='int16'),
        'satzenith': np.array([0, 9000], dtype='int16'),
        'azimuthdiff': np.array([0, 18000], dtype='int16'),
        'sunazimuth': np.array([-18000, 18000], dtype='int16'),
        'satazimuth': np.array([-18000, 18000], dtype='int16'),
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
    enc = {}
    if id_tag is not None:
        if id_tag.startswith('ch_tb'):
            # IR channel
            enc = {'dtype': 'int16',
                   'scale_factor': 0.01,
                   '_FillValue': -32767,
                   'zlib': True,
                   'complevel': 4,
                   'add_offset': 273.15}
        elif id_tag.startswith('ch_r'):
            # Refl channel
            enc = {'dtype': 'int16',
                   'scale_factor': 0.01,
                   'zlib': True,
                   'complevel': 4,
                   '_FillValue': -32767,
                   'add_offset': 0.0}
        elif id_tag in PPS_ANGLE_TAGS:
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
    if name in ['lon', 'lat']:
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
    if not enc:
        raise ValueError('Unsupported band: {}'.format(name))
    return name, enc


def remove_attributes(scene, band, remove):
    """Remove attributes from band."""
    for attr in remove:
        scene[band].attrs.pop(attr, None)


def rename_latitude_longitude(scene):
    """Rename latitude longitude to lat lon."""
    lat_name_satpy = 'latitude'
    lon_name_satpy = 'longitude'
    for alt_latname in ['lat_pixels', 'm_latitude', 'i_latitude']:
        if alt_latname in scene and 'latitude' not in scene:
            lat_name_satpy = alt_latname
    for alt_lonname in ['lon_pixels', 'm_longitude', 'i_longitude']:
        if alt_lonname in scene and 'longitude' not in scene:
            lon_name_satpy = alt_lonname
    # scene[lat_name_satpy].rename('lat')
    # scene[lon_name_satpy].rename('lon')
    scene[lat_name_satpy].attrs['name'] = 'lat'
    scene[lon_name_satpy].attrs['name'] = 'lon'
    scene['lat'] = scene[lat_name_satpy]
    scene['lon'] = scene[lon_name_satpy]
    del scene[lat_name_satpy]
    del scene[lon_name_satpy]

    # Update attributes
    scene['lat'].attrs['long_name'] = 'latitude coordinate'
    scene['lon'].attrs['long_name'] = 'longitude coordinate'
    scene['lat'].attrs['name'] = 'lat'
    scene['lon'].attrs['name'] = 'lon'
    scene['lon'].attrs['name'] = 'lon'
    scene['lon'].attrs['name'] = 'lon'
    scene['lon'].attrs['valid_range'] = np.array([-18000, 18000], dtype='float32')
    scene['lat'].attrs['valid_range'] = np.array([-9000, 90000], dtype='float32')
    for attr in ['valid_min', 'valid_max', 'coordinates',
                 'resolution', 'calibration', 'polarization', 'level',
                 'modifiers', '_satpy_id']:
        scene['lat'].attrs.pop(attr, None)
        scene['lon'].attrs.pop(attr, None)
    for coord_name in ['acq_time', 'm_latitude', 'i_latitude', 'm_latitude', 'i_latitude', 'latitude', 'longitude']:
        try:
            del scene['lat'].coords[coord_name]
            del scene['lon'].coords[coord_name]
        except KeyError:
            pass


def set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch):
    """Add some default values for band attributes."""
    nimg = 0  # name of first dataset is image0
    # Set some header attributes:
    scene.attrs['history'] = "Created by level1c4pps."
    if 'platform_name' in irch.attrs and 'platform' not in scene.attrs:
        scene.attrs['platform'] = irch.attrs['platform_name']
    if 'platform' in irch.attrs and 'platform' not in scene.attrs:
        scene.attrs['platform'] = irch.attrs['platform']
    if 'sensor' in irch.attrs:  # prefer channel sensor (often one)
        sensor_name = irch.attrs['sensor']
    elif 'sensor' in scene.attrs:  # might be a list
        if isinstance(scene.attrs['sensor'], list):
            sensor_name = scene.attrs['sensor'][0]
        else:
            sensor_name = scene.attrs['sensor']
    scene.attrs['sensor'] = sensor_name.upper()
    scene.attrs['instrument'] = sensor_name.upper()
    nowutc = datetime.utcnow()
    scene.attrs['date_created'] = nowutc.strftime("%Y-%m-%dT%H:%M:%SZ")
    scene.attrs['version_level1c4pps_satpy'] = satpy.__version__
    scene.attrs['version_level1c4pps'] = level1c4pps.__version__
    # bands
    for band in BANDNAMES:
        if band not in scene:
            continue
        idtag = PPS_TAGNAMES.get(band, None)
        if idtag is not None:
            scene[band].attrs['id_tag'] = idtag
        scene[band].attrs['description'] = sensor_name.upper() + ' ' + str(band).upper()
        if 'sun_earth_distance_correction_factor' not in scene[band].attrs.keys():
            scene[band].attrs['sun_earth_distance_correction_factor'] = 1.0
            scene[band].attrs['sun_earth_distance_correction_applied'] = 'False'
        else:
            # Assume factor applied if available as attribute.
            scene[band].attrs['sun_earth_distance_correction_applied'] = 'True'
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
        scene[band].attrs['name'] = "image{:d}".format(nimg)
        scene[band].attrs['coordinates'] = 'lon lat'
        if band in REFL_BANDS:
            scene[band].attrs['valid_range'] = np.array([0, 20000], dtype='int16')
            scene[band].attrs['units'] = '%'  # Needed by AVHRR
        else:
            scene[band].attrs['valid_range'] = np.array([-273.15 * 100, 300 * 100], dtype='int16')
            scene[band].attrs['units'] = 'K'  # Needed by AVHRR

        # Add time coordinate. To make cfwriter aware that we want 3D data.
        scene[band].coords['time'] = irch.attrs['start_time']

        # Remove some attributes and coordinates
        for attr in ['area', 'valid_min', 'valid_max']:
            scene[band].attrs.pop(attr, None)
        for coord_name in ['acq_time', 'latitude', 'longitude']:
            try:
                del scene[band].coords[coord_name]
            except KeyError:
                pass

        nimg += 1
    return nimg


def update_angle_attributes(scene, band):
    """Set and delete angle attributes."""
    for angle in PPS_ANGLE_TAGS:
        if angle not in scene and angle in ['sunazimuth', 'satazimuth']:
            # azimuth angles not always there
            continue
        scene[angle].attrs = {}
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
        for attr in ['area', 'valid_min', 'valid_max']:
            scene[angle].attrs.pop(attr, None)
            try:
                del scene[angle].encoding['coordinates']
            except (AttributeError, KeyError):
                pass
        # delete some coords
        for coord_name in ['acq_time', 'latitude', 'longitude']:
            try:
                del scene[angle].coords[coord_name]
            except KeyError:
                pass


def apply_sunz_correction(scene, REFL_BANDS):
    """Apply sun zenith angle correciton to visual channels."""
    sza = scene['sunzenith']
    mu0 = np.cos(np.radians(sza))
    scaler = 24.35 / (2 * mu0 + np.sqrt(498.5225 * mu0 * mu0 + 1))
    for band in REFL_BANDS:
        if band not in scene:
            continue
        if scene[band].attrs['sun_zenith_angle_correction_applied'] == 'False':
            scene[band].values = scene[band].values * scaler
            scene[band].attrs['sun_zenith_angle_correction_applied'] = 'True'


def platform_name_to_use_in_filename(platform_name):
    """Get platform name for PPS filenames from platfrom attribute."""
    new_name = platform_name.lower().replace('-', '').replace('aqua', '2').replace('terra', '1').replace("suomi", "")
    return new_name


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
            platform_name_to_use_in_filename(platform_name),
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
