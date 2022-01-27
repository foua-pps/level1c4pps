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
import satpy
from satpy.scene import Scene
from level1c4pps import (get_encoding, compose_filename,
                         set_header_and_band_attrs_defaults,
                         rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         convert_angles)
import pyspectral  # testing that pyspectral is available # noqa: F401
import logging
from satpy.utils import debug_on

from distutils.version import LooseVersion
if LooseVersion(satpy.__version__) < LooseVersion('0.22.1'):
    raise ImportError("'slstr2pps' requires satpy 0.22.1 or greater")

debug_on()

# Example:
# S3A_SL_1_RBT____20200911T163815_20200911T164115_20200912T215136_0179_062_340_0360_LN2_O_NT_004.SEN3#

logger = logging.getLogger('slstr2pps')

# slstr2pps.py -o temp/ S3A_SL_1_RBT____20200911T163815_*.SEN3/*n.nc
# ['F1_in', 'F2_in', 'S1_an', 'S2_an', 'S3_an', 'S4_an', 'S4_bn', 'S5_an', 'S5_bn',
# 'S6_an', 'S6_bn', 'S7_in', 'S8_in', 'S9_in',
# 'bayes_an', 'bayes_bn', 'bayes_in', 'cloud_an', 'cloud_bn', 'cloud_in',
# 'confidence_an', 'confidence_bn', 'confidence_in',
# 'latitude_an', 'latitude_ao', 'latitude_bn', 'latitude_bo', 'latitude_in', 'latitude_io',
# 'longitude_an', 'longitude_ao', 'longitude_bn', 'longitude_bo', 'longitude_in', 'longitude_io',
# 'pointing_an', 'pointing_bn', 'pointing_in',
# 'satellite_azimuth_angle_n', 'satellite_azimuth_angle_o', 'satellite_zenith_angle_n', 'satellite_zenith_angle_o',
# 'solar_azimuth_angle_n', 'solar_azimuth_angle_o', 'solar_zenith_angle_n', 'solar_zenith_angle_o']


BANDNAMES = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'F1', 'F2']

REFL_BANDS = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']

ANGLE_NAMES = ['satellite_azimuth_angle', 'satellite_zenith_angle',
               'solar_azimuth_angle', 'solar_zenith_angle']

PPS_TAGNAMES = {'S2': 'ch_r06',  # or S1
                'S3': 'ch_r09',
                'S4': 'ch_r13',
                'S5': 'ch_r16',
                'S7': 'ch_tb37',
                'S8': 'ch_tb11',
                'S9': 'ch_tb12',
                # Not yet in pps:
                'S6': 'ch_r21',
                'S1': 'ch_rxx',
                'F1': 'ch_tbxx',
                'F2': 'ch_tbxx'}

BANDNAMES_PPS = ['S2', 'S3', 'S4', 'S5', 'S7', 'S8', 'S9']
BANDNAMES_DEFAULT = ['S2', 'S3', 'S4', 'S5',  'S6', 'S7', 'S8', 'S9']


def get_encoding_slstr(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene, orbit_n=0):
    """Set and delete some attributes."""
    irch = scene['S8']
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "slstr2pps.py"
    return nimg


def process_one_scene(scene_files, out_path, engine='h5netcdf',
                      all_channels=False, pps_channels=False, orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(
        reader='slstr_l1b',
        filenames=scene_files)

    MY_BANDNAMES = BANDNAMES_DEFAULT
    if all_channels:
        MY_BANDNAMES = BANDNAMES
    if pps_channels:
        MY_BANDNAMES = BANDNAMES_PPS

    scn_.load(MY_BANDNAMES + ['latitude', 'longitude'] + ANGLE_NAMES)

    # Everything should be on the same grid, to be saved as ppsleve1c
    scn_ = scn_.resample(resampler="native")

    # one ir channel
    irch = scn_['S8']

    # Set header and band attributes
    set_header_and_band_attrs(scn_, orbit_n=orbit_n)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, delete_azimuth=True)
    update_angle_attributes(scn_, irch)

    filename = compose_filename(scn_, out_path, instrument='slstr', band=irch)
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='slstr'),
                       engine=engine,
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_slstr(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
