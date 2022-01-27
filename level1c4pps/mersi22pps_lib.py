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

"""Functions to convert MERSI-2 level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import os
import time
from satpy.scene import Scene
from level1c4pps import (get_encoding, compose_filename,
                         set_header_and_band_attrs_defaults,
                         ANGLE_ATTRIBUTES, rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         convert_angles)
import pyspectral  # testing that pyspectral is available # noqa: F401
import logging

# Example:
# tf2019234102243.FY3D-X_MERSI_GEOQK_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_GEO1K_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_1000M_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_0250M_L1B.HDF
#

logger = logging.getLogger('mersi22pps')

PLATFORM_SHORTNAMES = {'FY3D': 'FY-3D'}
# BANDNAMES = ['%d' % (chn+1) for chn in range(25)]
BANDNAMES = ['3', '4', '5', '6', '20', '23', '24', '25']

REFL_BANDS = ['3', '4', '5', '6']

ANGLE_NAMES = ['satellite_zenith_angle', 'solar_zenith_angle',
               'satellite_azimuth_angle', 'solar_azimuth_angle']

PPS_TAGNAMES = {'3': 'ch_r06',
                '4': 'ch_r09',
                '5': 'ch_r13',
                '6': 'ch_r16',
                '20': 'ch_tb37',
                '22': 'ch_tb73',
                '23': 'ch_tb85',
                '24': 'ch_tb11',
                '25': 'ch_tb12'}

MERSI2_LEVEL1_FILE_PATTERN = 'tf{start_time:%Y%j%H%M%S}.{platform_shortname:4s}-X_MERSI_{dataset}_L1B.HDF'


def get_encoding_mersi2(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene, orbit_n=0):
    """Set and delete some attributes."""
    irch = scene['24']
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "mersi22pps.py"
    return nimg


def remove_broken_data(scene):
    """Set bad data to nodata."""
    import numpy as np
    for band in BANDNAMES:
        if band in REFL_BANDS:
            continue
        if band in scene:
            remove = np.where(scene[band].values < 1, np.nan, 0)  # 1K very cold
            scene[band].values = scene[band].values + remove


def process_one_scene(scene_files, out_path, engine='h5netcdf', orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(
        reader='mersi2_l1b',
        filenames=scene_files)

    scn_.load(BANDNAMES + ['latitude', 'longitude'] + ANGLE_NAMES, resolution=1000)

    # Remove bad data at first and last column
    remove_broken_data(scn_)

    # one ir channel
    irch = scn_['24']

    # Set header and band attributes
    set_header_and_band_attrs(scn_, orbit_n=orbit_n)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, delete_azimuth=True)
    update_angle_attributes(scn_, irch)
    for angle in ['sunzenith', 'satzenith', 'azimuthdiff']:
        scn_[angle].attrs['file_key'] = ANGLE_ATTRIBUTES['mersi2_file_key'][angle]

    filename = compose_filename(scn_, out_path, instrument='mersi2', band=irch)
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='mersi-2'),
                       engine=engine,
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_mersi2(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
    return filename
