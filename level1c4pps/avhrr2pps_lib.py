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

"""Functions to convert AVHRR AAPP or EPS l1b data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import os
import time
from satpy.scene import Scene
from level1c4pps import (get_encoding, compose_filename,
                         rename_latitude_longitude,
                         set_header_and_band_attrs_defaults,
                         update_angle_attributes, get_header_attrs,
                         convert_angles,
                         apply_sunz_correction)

import logging

# Example:
# hrpt_noaa19_20141210_1056_30086.l1b

logger = logging.getLogger('avhrr2pps')

BANDNAMES = ['1', '2', '3a', '3b', '4', '5']

ANGLE_NAMES_EPS = ['satellite_zenith_angle', 'solar_zenith_angle',
                   'satellite_azimuth_angle', 'solar_azimuth_angle']
ANGLE_NAMES_AAPP = ['sensor_zenith_angle', 'solar_zenith_angle',
                    'sun_sensor_azimuth_difference_angle']

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


def set_header_and_band_attrs(scene, orbit_n=0):
    """Set and delete some attributes."""
    irch = scene['4']
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "avhrr2pps.py"
    return nimg


def check_broken_data(scene):
    """Set bad data to nodata."""
    import numpy as np
    lat = scene['latitude']
    # If we have data in line 2 it is ok
    if (lat[1, :].values > 0).any():
        return
    lon = scene['longitude']
    part_of_lat_that_is_zero = np.sum(lat.values.ravel() == 0) / len(lat.values.ravel())
    part_of_lon_that_is_zero = np.sum(lon.values.ravel() == 0) / len(lon.values.ravel())
    if (part_of_lat_that_is_zero > 0.90 and part_of_lon_that_is_zero > 0.90):
        raise ValueError(
            'More than 90% of data have (lat, lon) equal to (0,0).\n'
            'Most likely this is an older format of hrpt not yet suported by satpy.\n'
            'Stopping. File will not be written.')


def process_one_scene(scene_files, out_path, engine='h5netcdf', orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    if 'AVHR_xxx' in scene_files[0]:
        avhrr_reader = 'avhrr_l1b_eps'
        angles = ANGLE_NAMES_EPS
    else:
        avhrr_reader = 'avhrr_l1b_aapp'
        angles = ANGLE_NAMES_AAPP
    scn_ = Scene(
        reader=avhrr_reader,
        filenames=scene_files)
    scn_.load(BANDNAMES + ['latitude', 'longitude'] + angles)
    # one ir channel
    irch = scn_['4']

    # Check if we have old hrpt format with data only every 20th line
    check_broken_data(scn_)

    # Set header and band attributes
    set_header_and_band_attrs(scn_, orbit_n=orbit_n)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, delete_azimuth=True)
    update_angle_attributes(scn_, irch)

    # Apply sunz correction
    apply_sunz_correction(scn_, REFL_BANDS)

    filename = compose_filename(scn_, out_path, instrument='avhrr', band=irch)
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='avhrr'),
                       engine=engine,
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_avhrr(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
    return filename
