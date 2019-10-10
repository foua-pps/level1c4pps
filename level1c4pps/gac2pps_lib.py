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
from datetime import datetime
from satpy.scene import Scene
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


def process_one_file(gac_file, out_path='.'):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    image_num = 0  # name of first dataset is image0
    # platform_shortname = p__.parse(
    #     os.path.basename(tslot_files[0]))['platform_shortname']
    # start_time = p__.parse(
    #     os.path.basename(tslot_files[0]))['start_time']
    # platform_name = PLATFORM_SHORTNAMES[platform_shortname]
    # #Load channel data for one scene and set some attributes
    # coefs = get_calibration_for_time(platform=platform_shortname,
    #                                  time=start_time)

    scn_ = Scene(
        reader='avhrr_l1b_gaclac',
        filenames=[gac_file])

    if 'avhrr-3' in scn_.attrs['sensor']:
        sensor = 'avhrr'
        scn_.load(BANDNAMES + ['latitude', 'longitude',
                               'sensor_zenith_angle', 'solar_zenith_angle',
                               'solar_azimuth_angle', 'sensor_azimuth_angle',
                               'sun_sensor_azimuth_difference_angle'])
    for band in BANDNAMES:
        try:
            idtag = PPS_TAGNAMES.get(band, band)
            scn_[band].attrs['id_tag'] = idtag
            scn_[band].attrs['description'] = 'AVHRR ' + str(band)
            scn_[band].attrs['sun_earth_distance_correction_applied'] = 'False'
            scn_[band].attrs['sun_earth_distance_correction_factor'] = 1.0
            scn_[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
            scn_[band].attrs['name'] = "image{:d}".format(image_num)
            scn_[band].attrs['coordinates'] = 'lon lat'
            del scn_[band].attrs['area']
            image_num += 1
        except KeyError:
            continue

    # Set some header attributes:
    scn_.attrs['instrument'] = sensor.upper()
    scn_.attrs['source'] = "gac2pps.py"
    nowutc = datetime.utcnow()
    scn_.attrs['date_created'] = nowutc.strftime("%Y-%m-%dT%H:%M:%SZ")

    # Find lat/lon data
    irch = scn_['4']
    scn_.attrs['platform'] = irch.attrs['platform_name']
    scn_.attrs['platform_name'] = irch.attrs['platform_name']
    scn_.attrs['orbit_number'] = '{:05d}'.format(irch.attrs['orbit_number'])
    scn_.attrs['orbit'] = scn_.attrs['orbit_number']
    # lons = lons.where(lons <= 360, -999.0)
    # lons = lons.where(lons >= -360, 999.0)
    # lats = lats.where(lats <= 90, -999.0)
    # lats = lats.where(lats >= -90, 999.0)

    scn_['lat'] = scn_['latitude']
    del scn_['latitude']
    scn_['lat'].attrs['long_name'] = 'latitude coordinate'
    del scn_['lat'].coords['acq_time']

    scn_['lon'] = scn_['longitude']
    del scn_['longitude']
    scn_['lon'].attrs['long_name'] = 'longitude coordinate'
    del scn_['lon'].coords['acq_time']

    angle_names = []
    scn_['sunzenith'] = scn_['solar_zenith_angle']
    del scn_['solar_zenith_angle']
    scn_['sunzenith'].attrs['id_tag'] = 'sunzenith'
    scn_['sunzenith'].attrs['long_name'] = 'sun zenith angle'
    scn_['sunzenith'].attrs['valid_range'] = [0, 18000]
    scn_['sunzenith'].attrs['name'] = "image{:d}".format(image_num)
    angle_names.append("image{:d}".format(image_num))
    scn_['sunzenith'].attrs['coordinates'] = 'lon lat'
    del scn_['sunzenith'].attrs['area']
    scn_['sunzenith'].coords['time'] = irch.attrs['start_time']
    del scn_['sunzenith'].coords['acq_time']
    image_num += 1

    # satzenith
    scn_['satzenith'] = scn_['sensor_zenith_angle']
    del scn_['sensor_zenith_angle']
    scn_['satzenith'].attrs['id_tag'] = 'satzenith'
    scn_['satzenith'].attrs['long_name'] = 'satellite zenith angle'
    scn_['satzenith'].attrs['valid_range'] = [0, 9000]
    scn_['satzenith'].attrs['name'] = "image{:d}".format(image_num)
    angle_names.append("image{:d}".format(image_num))
    scn_['satzenith'].attrs['coordinates'] = 'lon lat'
    del scn_['satzenith'].attrs['area']
    scn_['satzenith'].coords['time'] = irch.attrs['start_time']
    del scn_['satzenith'].coords['acq_time']
    image_num += 1

    # azidiff
    scn_['azimuthdiff'] = abs(scn_['sun_sensor_azimuth_difference_angle'])
    scn_['azimuthdiff'].attrs = scn_['sun_sensor_azimuth_difference_angle'].attrs
    del scn_['sun_sensor_azimuth_difference_angle']
    scn_['azimuthdiff'].attrs['id_tag'] = 'azimuthdiff'
    # scn_['azimuthdiff'].attrs['standard_name'] = (
    #     'angle_of_rotation_from_solar_azimuth_to_platform_azimuth')
    scn_['azimuthdiff'].attrs['long_name'] = 'absolute azimuth difference angle'
    scn_['azimuthdiff'].attrs['valid_range'] = [0, 18000]
    scn_['azimuthdiff'].attrs['name'] = "image{:d}".format(image_num)
    angle_names.append("image{:d}".format(image_num))
    scn_['azimuthdiff'].attrs['coordinates'] = 'lon lat'
    del scn_['azimuthdiff'].attrs['area']
    scn_['azimuthdiff'].coords['time'] = irch.attrs['start_time']
    del scn_['azimuthdiff'].coords['acq_time']
    image_num += 1

    # satazimuth
    scn_['sunazimuth'] = scn_['solar_azimuth_angle']
    del scn_['solar_azimuth_angle']
    scn_['sunazimuth'].attrs['id_tag'] = 'sunazimuth'
    scn_['sunazimuth'].attrs['long_name'] = 'sun azimuth angle'
    scn_['sunazimuth'].attrs['valid_range'] = [0, 18000]
    scn_['sunazimuth'].attrs['name'] = "image{:d}".format(image_num)
    angle_names.append("image{:d}".format(image_num))
    scn_['sunazimuth'].attrs['coordinates'] = 'lon lat'
    del scn_['sunazimuth'].attrs['area']
    scn_['sunazimuth'].coords['time'] = irch.attrs['start_time']
    del scn_['sunazimuth'].coords['acq_time']
    image_num += 1

    # satazimuth
    scn_['satazimuth'] = scn_['sensor_azimuth_angle']
    del scn_['sensor_azimuth_angle']
    scn_['satazimuth'].attrs['id_tag'] = 'satazimuth'
    scn_['satazimuth'].attrs['long_name'] = 'satellite azimuth angle'
    scn_['satazimuth'].attrs['valid_range'] = [0, 9000]
    scn_['satazimuth'].attrs['name'] = "image{:d}".format(image_num)
    angle_names.append("image{:d}".format(image_num))
    scn_['satazimuth'].attrs['coordinates'] = 'lon lat'
    del scn_['satazimuth'].attrs['area']
    scn_['satazimuth'].coords['time'] = irch.attrs['start_time']
    del scn_['satazimuth'].coords['acq_time']
    image_num += 1

    # Get filename
    start_time = irch.attrs['start_time']
    end_time = irch.attrs['end_time']
    platform_name = irch.attrs['platform_name']
    orbit_number = int(scn_.attrs['orbit_number'])
    filename = os.path.join(
        out_path,
        "S_NWC_avhrr_{:s}_{:05d}_{:s}Z_{:s}Z.nc".format(
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
        del scn_[band].coords['acq_time']

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
        save_info[name] = {'dtype': 'int16',
                           'scale_factor': 0.01,
                           'zlib': True,
                           'complevel': 4,
                           '_FillValue': -32767,
                           'add_offset': 0.0}
    for name in ['lon', 'lat']:
        save_info[name] = {'dtype': 'float32',    'zlib': True,
                           'complevel': 4, '_FillValue': -999.0}
    header_attrs = scn_.attrs.copy()
    header_attrs['start_time'] = time.strftime(
        "%Y-%m-%d %H:%M:%S",
        irch.attrs['start_time'].timetuple())
    header_attrs['end_time'] = time.strftime(
        "%Y-%m-%d %H:%M:%S",
        irch.attrs['end_time'].timetuple())
    header_attrs['sensor'] = sensor.lower()

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
