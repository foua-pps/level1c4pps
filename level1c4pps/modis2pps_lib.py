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

"""Functions to convert MODIS level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import os
import time
from datetime import datetime
from satpy.scene import Scene
from level1c4pps import (get_encoding, compose_filename,
                         rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         convert_angles,
                         apply_sunz_correction)

import logging

# Example:
# MYD03.A2014272.0215.006.2014274124038.hdf
# MYD021KM.A2014272.0215.006.2014274142955.hdf


logger = logging.getLogger('modis2pps')

BANDNAMES = ['1', '2', '6', '20', '26', '28', '29', '31', '32']

REFL_BANDS = ['1', '2', '6', '26']

ANGLE_NAMES = ['satellite_zenith_angle', 'solar_zenith_angle',
               'satellite_azimuth_angle', 'solar_azimuth_angle']

PPS_TAGNAMES = {'1':  'ch_r06',
                '2':  'ch_r09',
                '26': 'ch_r13',
                '6':  'ch_r16',
                '20': 'ch_tb37',
                '28': 'ch_tb73',
                '29': 'ch_tb85',
                '31': 'ch_tb11',
                '32': 'ch_tb12'}


def get_encoding_modis(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene):
    """Set and delete some attributes."""
    nimg = 0  # name of first dataset is image0
    # Set some header attributes:
    irch = scene['31']
    scene.attrs['platform'] = irch.attrs['platform_name']
    sensor_name = [x for x in scene.attrs['sensor']][0]
    scene.attrs['instrument'] = sensor_name.upper()
    scene.attrs['source'] = "modis2pps.py"
    nowutc = datetime.utcnow()
    scene.attrs['date_created'] = nowutc.strftime("%Y-%m-%dT%H:%M:%SZ")
    scene.attrs['orbit_number'] = 00000
    # bands
    for band in BANDNAMES:
        if band not in scene:
            continue
        idtag = PPS_TAGNAMES.get(band, None)
        if not idtag:
            continue
        scene[band].attrs['id_tag'] = idtag
        scene[band].attrs['description'] = 'MODIS ' + str(band)
        scene[band].attrs['sun_earth_distance_correction_applied'] = 'False'
        scene[band].attrs['sun_earth_distance_correction_factor'] = 1.0
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
        scene[band].attrs['name'] = "image{:d}".format(nimg)
        scene[band].attrs['coordinates'] = 'lon lat'
        # Add time coordinate. To make cfwriter aware that we want 3D data.
        scene[band].coords['time'] = irch.attrs['start_time']
        try:
            del scene[band].attrs['area']
        except KeyError:
            pass
        nimg += 1
    return nimg


def process_one_scene(scene_files, out_path):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(
        reader='modis_l1b',
        filenames=scene_files)

    scn_.load(BANDNAMES + ['latitude', 'longitude'] + ANGLE_NAMES, resolution=1000)
    # one ir channel
    irch = scn_['31']

    # Remove som buggy channels:
    if 'aqua' in irch.attrs['platform_name'].lower():
        del(scn_['6'])  # Skip 1.6 for Aqua
    elif 'terra' in irch.attrs['platform_name'].lower():
        del(scn_['29'])  # Skip 8.6 for Terra

    # Set header and band attributes
    set_header_and_band_attrs(scn_)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, delete_azimuth=True)
    update_angle_attributes(scn_, irch)

    # Apply sunz correction
    apply_sunz_correction(scn_, REFL_BANDS)

    filename = compose_filename(scn_, out_path, instrument='modis', band=irch)
    filename = filename.replace('aqua', '2')
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='modis'),
                       engine='netcdf4',
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_modis(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
