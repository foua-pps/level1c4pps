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

"""Functions to convert FIDUCEO level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import os
import time
from satpy.scene import Scene
from level1c4pps import (get_encoding, compose_filename,
                         set_header_and_band_attrs_defaults,
                         rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         convert_angles,
                         apply_sunz_correction)

import logging

# FIDUCEO_FCDR_L1C_AVHRR_N18ALL_20080718161223_20080718175428_EASY_v1.00_fv2.0.0.nc

logger = logging.getLogger('fiduceo2pps')

BANDNAMES = ['1', '2', '3', '3a', '3b', '4', '5']

REFL_BANDS = ['1', '2', '3a']

PPS_TAGNAMES = {"1": "ch_r06",
                "2": "ch_r09",
                "3a": "ch_r16",
                "3": "ch_tb37",
                "3b": "ch_tb37",
                "4": "ch_tb11",
                "5": "ch_tb12"}

ANGLE_NAMES = ['sensor_zenith_angle', 'solar_zenith_angle',
               'sun_sensor_azimuth_difference_angle']

MOVE_TO_HEADER_ATTRS = ['source',
                        'Ch3a_only',
                        'file_type',
                        'Ch3a_Ch3b_split_file',
                        'licence',
                        'naming_authority',
                        'software_version',
                        'UUID',
                        'writer_version',
                        'template_key',
                        'Conventions',
                        'institution',
                        'comment',
                        'references',
                        'history',
                        'title']


def get_encoding_fiduceo(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene):
    """Set and delete some attributes."""
    irch = scene['4']
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch)
    scene.attrs['source'] = "fiduceo2pps.py"
    scene.attrs['orbit_number'] = 99999
    # Move some attributes to header:
    for attr in MOVE_TO_HEADER_ATTRS:
        scene.attrs[attr] = irch.attrs[attr]
    for band in BANDNAMES:
        if band not in scene:
            continue
        if band in REFL_BANDS:
            # For FIDUCEO data sun_earth_distance_correction is applied always!
            scene[band].attrs['sun_earth_distance_correction_applied'] = 'True'
    return nimg


def remove_some_attrs_and_none(scene):
    """Remove some attributes and replace some None with none."""
    for dataset in scene:
        print(dataset.name)
        for coord in dataset.coords:
            for attr in MOVE_TO_HEADER_ATTRS:
                dataset.coords[coord].attrs.pop(attr, None)
            for attr in dataset.coords[coord].attrs:
                if dataset.coords[coord].attrs[attr] is None and attr not in ['_FillValue']:
                    dataset.coords[coord].attrs[attr] = 'none'
        for attr in MOVE_TO_HEADER_ATTRS:
            dataset.attrs.pop(attr, None)
        for attr in dataset.attrs:
            if dataset.attrs[attr] is None and attr not in ['_FillValue']:
                dataset.attrs[attr] = 'none'


def process_one_scene(scene_files, out_path):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(
        reader='fiduceo_avhrr_fcdr_nc',
        filenames=scene_files)
    scn_.load(BANDNAMES + ['latitude', 'longitude'] + ANGLE_NAMES)
    # one ir channel
    irch = scn_['4']

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Set header and band attributes
    set_header_and_band_attrs(scn_)

    # Convert angles to PPS
    convert_angles(scn_, delete_azimuth=False)
    update_angle_attributes(scn_, irch)

    # Apply sunz correction
    apply_sunz_correction(scn_, REFL_BANDS)

    remove_some_attrs_and_none(scn_)
    filename = compose_filename(scn_, out_path, instrument='avhrr', band=irch)
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='avhrr'),
                       engine='netcdf4',
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_fiduceo(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
