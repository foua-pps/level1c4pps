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

"""Functions to convert AVHRR HRPT or EPS l1b data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import os
import time
from datetime import datetime
from satpy.scene import Scene
from level1c4pps import (get_encoding, compose_filename,
                         rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         SATPY_ANGLE_NAMES, convert_angles,
                         apply_sunz_correction)

import logging

# Example:
# hrpt_noaa19_20141210_1056_30086.l1b

logger = logging.getLogger('avhrr2pps')

BANDNAMES = ['1', '2', '3a', '3b', '4', '5']

REFL_BANDS = ['1', '2', '3a']

PPS_TAGNAMES = {'1':  'ch_r06',
                '2':  'ch_r09',
                '3a':  'ch_r16',
                '3b': 'ch_tb37',
                '4': 'ch_tb11',
                '5': 'ch_tb12'}


def get_encoding_avhrr(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene):
    """Set and delete some attributes."""
    nimg = 0  # name of first dataset is image0
    # Set some header attributes:
    irch = scene['4']
    scene.attrs['platform'] = irch.attrs['platform_name']
    sensor_name = [x for x in scene.attrs['sensor']][0]
    scene.attrs['instrument'] = sensor_name.upper()
    scene.attrs['source'] = "avhrr2pps.py"
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
        scene[band].attrs['description'] = 'AVHRR ' + str(band)
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
    if 'AVHR_xxx' in scene_files[0]:
        avhrr_reader = 'avhrr_l1b_eps'
    else:
        avhrr_reader = 'aapp_l1b'   
    scn_ = Scene(
        reader=avhrr_reader,
        filenames=scene_files)
    #import pdb;pdb.set_trace()
    scn_.load(BANDNAMES + ['latitude', 'longitude',
                           'satellite_zenith_angle', 'solar_zenith_angle',
                           'satellite_azimuth_angle', 'solar_azimuth_angle'])
    # one ir channel
    irch = scn_['4']

    # Set header and band attributes
    set_header_and_band_attrs(scn_)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, SATPY_ANGLE_NAMES)
    update_angle_attributes(scn_, irch)

    # Apply sunz correction
    apply_sunz_correction(scn_, REFL_BANDS)

    filename = compose_filename(scn_, out_path, instrument='avhrr', band=irch)
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='avhrr'),
                       engine='netcdf4',
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_avhrr(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
