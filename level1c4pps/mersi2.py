#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c21529.ad.smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Functions to convert MERSI-2 level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file
"""

import os
import numpy as np
import xarray as xr
import dask.array as da
from glob import glob
import time
from datetime import datetime
from satpy.scene import Scene
from trollsift.parser import globify, Parser
from pyorbital.astronomy import get_alt_az, sun_zenith_angle
from pyorbital.orbital import get_observer_look
from level1c4pps.calibration_coefs import get_calibration_for_time, CALIB_MODE
from level1c4pps import make_azidiff_angle
import pyresample
import logging

# Example:
# tf2019234102243.FY3D-X_MERSI_GEOQK_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_GEO1K_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_1000M_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_0250M_L1B.HDF
#

logger = logging.getLogger('mersi22pps')

PLATFORM_SHORTNAMES = {'FY3D': 'FY-3D'}
#BANDNAMES = ['%d' % (chn+1) for chn in range(25)]
BANDNAMES = ['5', '6', '12', '15', '20', '23', '24', '25']

PPS_TAGNAMES = {'12': 'ch_r06',
                '15': 'ch_r09',
                '5': 'ch_r13',
                '6': 'ch_r16',
                '20': 'ch_tb37',
                '23': 'ch_tb85',
                '24': 'ch_tb11',
                '25': 'ch_tb12'}


MERSI2_LEVEL1_FILE_PATTERN = 'tf{start_time:%Y%j%H%M%S}.{platform_shortname:4s}-X_MERSI_{dataset}_L1B.HDF'

p__ = Parser(MERSI2_LEVEL1_FILE_PATTERN)


def process_one_scene(scene_files, out_path):
    """ Make level 1c files in PPS-format """
    tic = time.time()

    image_num = 0  # name of first dataset is image0

    platform_shortname = p__.parse(
        os.path.basename(scene_files[0]))['platform_shortname']
    start_time = p__.parse(
        os.path.basename(scene_files[0]))['start_time']
    platform_name = PLATFORM_SHORTNAMES[platform_shortname]

    scn_ = Scene(
        reader='mersi2_l1b',
        filenames=scene_files)

    scn_.attrs['platform_name'] = platform_name
    scn_.load(BANDNAMES + ['latitude', 'longitude',
                           'satellite_zenith_angle', 'solar_zenith_angle',
                           'satellite_azimuth_angle', 'solar_azimuth_angle'], resolution=1000)
    nimg = 0
    for band in BANDNAMES:
        idtag = PPS_TAGNAMES.get(band, None)
        if not idtag:
            continue

        scn_[band].attrs['id_tag'] = idtag
        scn_[band].attrs['description'] = 'MERSI-2 band ' + str(band)
        scn_[band].attrs['sun_earth_distance_correction_applied'] = 'False'
        scn_[band].attrs['sun_earth_distance_correction_factor'] = 1.0
        scn_[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
        scn_[band].attrs['name'] = "image{:d}".format(nimg)
        scn_[band].attrs['coordinates'] = 'lon lat'
        nimg += 1

    # Set som header attributes:
    scn_.attrs['platform'] = platform_name
    sensor_name = [x for x in scn_.attrs['sensor']][0]
    scn_.attrs['instrument'] = sensor_name.upper()
    scn_.attrs['source'] = "mersi22pps.py"
    # Perhaps one can get the orbit number from the hdf file?
    # FIXME!
    scn_.attrs['orbit_number'] = "99999"
    #scn_.attrs['orbit'] = scn_.attrs['orbit_number']
    irch = scn_['24']

    nowutc = datetime.utcnow()
    scn_.attrs['date_created'] = nowutc.strftime("%Y-%m-%dT%H:%M:%SZ")

    scn_['lat'] = scn_['latitude'].copy()
    del scn_['latitude']
    scn_['lat'].attrs['long_name'] = 'latitude coordinate'

    scn_['lon'] = scn_['longitude'].copy()
    del scn_['longitude']
    scn_['lon'].attrs['long_name'] = 'longitude coordinate'

    angle_names = []
    scn_['sunzenith'] = scn_['solar_zenith_angle']
    del scn_['solar_zenith_angle']
    scn_['sunzenith'].attrs['id_tag'] = 'sunzenith'
    scn_['sunzenith'].attrs['file_key'] = 'Geolocation/SolarZenithAngle'
    scn_['sunzenith'].attrs['long_name'] = 'solar zenith angle'
    scn_['sunzenith'].attrs['standard_name'] = 'solar_zenith_angle'
    scn_['sunzenith'].attrs['valid_range'] = [0, 18000]
    scn_['sunzenith'].attrs['name'] = "image{:d}".format(nimg)
    angle_names.append("image{:d}".format(image_num))
    scn_['sunzenith'].attrs['coordinates'] = 'lon lat'
    del scn_['sunzenith'].attrs['area']
    scn_['sunzenith'].coords['time'] = irch.attrs['start_time']
    nimg += 1

    # satzenith
    scn_['satzenith'] = scn_['satellite_zenith_angle']
    del scn_['satellite_zenith_angle']
    scn_['satzenith'].attrs['id_tag'] = 'satzenith'
    scn_['satzenith'].attrs['file_key'] = 'Geolocation/SensorZenithAngle'
    scn_['satzenith'].attrs['long_name'] = 'satellite zenith angle'
    scn_['satzenith'].attrs['stadard_name'] = 'platform_zenith_angle'
    scn_['satzenith'].attrs['valid_range'] = [0, 9000]
    scn_['satzenith'].attrs['name'] = "image{:d}".format(nimg)
    angle_names.append("image{:d}".format(nimg))
    scn_['satzenith'].attrs['coordinates'] = 'lon lat'
    del scn_['satzenith'].attrs['area']
    scn_['satzenith'].coords['time'] = irch.attrs['start_time']
    nimg += 1

    # azidiff
    scn_['azimuthdiff'] = make_azidiff_angle(scn_['satellite_azimuth_angle'], scn_['solar_azimuth_angle'])
    scn_['azimuthdiff'].attrs = scn_['solar_azimuth_angle'].attrs.copy()
    del scn_['solar_azimuth_angle']
    del scn_['satellite_azimuth_angle']
    scn_['azimuthdiff'].attrs['id_tag'] = 'azimuthdiff'
    scn_['azimuthdiff'].attrs['file_key'] = 'Geolocation/SensorSolarAzimuthDifference'
    scn_['azimuthdiff'].attrs['standard_name'] = 'absolute_sensor_solar_azimuth_difference_angle'
    scn_['azimuthdiff'].attrs['long_name'] = 'absolute sensor-solar azimuth difference angle'
    scn_['azimuthdiff'].attrs['Description'] = ('Sensor Solar azimuth diference angle at' +
                                                ' the geolocated beam position center')
    scn_['azimuthdiff'].attrs['valid_range'] = [0, 18000]
    scn_['azimuthdiff'].attrs['name'] = "image{:d}".format(nimg)
    angle_names.append("image{:d}".format(nimg))
    scn_['azimuthdiff'].attrs['coordinates'] = 'lon lat'
    del scn_['azimuthdiff'].attrs['area']
    scn_['azimuthdiff'].coords['time'] = irch.attrs['start_time']
    nimg += 1

    # Get filename
    start_time = irch.attrs['start_time']
    end_time = irch.attrs['end_time']
    platform_name = irch.attrs['platform_name']
    orbit_number = int(scn_.attrs['orbit_number'])

    filename = os.path.join(
        out_path,
        "S_NWC_mersi2_{:s}_{:05d}_{:s}Z_{:s}Z.nc".format(
            platform_name.lower().replace('-', ''),
            orbit_number,
            start_time.strftime('%Y%m%dT%H%M%S%f')[:-5],
            end_time.strftime('%Y%m%dT%H%M%S%f')[:-5]))

    for dataset in scn_.keys():
        if hasattr(scn_[dataset], 'attrs'):
            if hasattr(scn_[dataset].attrs, 'modifiers'):
                scn_[dataset].attrs['modifiers'] = 0.0

    # Encoding for channels
    save_info = {}
    for band in BANDNAMES:
        idtag = PPS_TAGNAMES[band]
        try:
            name = scn_[band].attrs['name']
        except KeyError:
            logger.debug("No band named %s", band)
            continue
        # Add time coordinate. To make cfwriter aware that we want 3D data.
        scn_[band].coords['time'] = irch.attrs['start_time']

        print(name)

        if 'tb' in idtag:
            save_info[name] = {'dtype': 'int16',
                               'scale_factor': 0.01,
                               '_FillValue': -32767,
                               'zlib': True,
                               'complevel': 4,
                               'add_offset': 273.15}
        else:
            save_info[name] = {'dtype': 'int16',
                               'scale_factor': 0.01,
                               'zlib': True,
                               'complevel': 4,
                               '_FillValue': -32767,
                               'add_offset': 0.0}

    # Encoding for angles and lat/lon
    for name in angle_names:
        print(name)
        save_info[name] = {'dtype': 'int16',
                           'scale_factor': 0.01,
                           'zlib': True,
                           'complevel': 4,
                           '_FillValue': -32767,
                           'add_offset': 0.0}
    for name in ['lon', 'lat']:
        print(name)
        save_info[name] = {'dtype': 'float32',    'zlib': True,
                           'complevel': 4, '_FillValue': -999.0}
    header_attrs = scn_.attrs.copy()
    header_attrs['start_time'] = time.strftime(
        "%Y-%m-%d %H:%M:%S",
        irch.attrs['start_time'].timetuple())
    header_attrs['end_time'] = time.strftime(
        "%Y-%m-%d %H:%M:%S",
        irch.attrs['end_time'].timetuple())
    header_attrs['sensor'] = sensor_name.lower()

    for band in BANDNAMES:
        idtag = PPS_TAGNAMES[band]
        try:
            to_pop = []
            for attr in scn_[band].attrs.keys():
                if hasattr(scn_[band].attrs[attr], 'keys'):
                    print("found dict", attr)
                    to_pop.append(attr)
            for attr in to_pop:
                attr_dict = scn_[band].attrs[attr]
                scn_[band].attrs.pop(attr)
                for key in attr_dict.keys():
                    scn_[band].attrs[attr+str(key)] = attr_dict[key]
        except KeyError:
            continue

    scn_.save_datasets(writer='cf', filename=filename,
                       header_attrs=header_attrs, engine='netcdf4',
                       encoding=save_info)
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
