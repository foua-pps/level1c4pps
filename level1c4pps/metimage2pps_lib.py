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

"""Functions to convert EPS-SG MetImage level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import numpy as np
import time
import satpy
from satpy.scene import Scene
from level1c4pps import (apply_sunz_correction,
                         compose_filename,
                         rename_latitude_longitude,
                         get_refl_bands,
                         get_band_names,
                         save_data,
                         check_file_exists,
                         log_time,
                         update_angle_attributes, get_header_attrs,
                         set_header_and_band_attrs_defaults,
                         convert_angles,
                         adjust_lons_to_valid_range)

import logging
from packaging.version import Version


# from satpy.utils import debug_on
# debug_on()

# Example:
# W_xx-eumetsat-darmstadt,SAT,SGA1-VII-1B-RAD_C_EUMT_20191001043852_G_D_20070912095840_20070912100343_T_N____.nc'
# W_XX-EUMETSAT-Darmstadt,SAT,SGA1-VII-1B-RAD_C_EUMT_20260211124622_G_V_20260211115100_20260211115158_C_N_T__.nc


logger = logging.getLogger('metimage2pps')

if Version(satpy.__version__) <= Version('0.59.0'):
    logger.info("For METimage use satpy >= 0.60.")
    raise ValueError("For METimage satpy >= 0.60 is needed")
else:
    logger.info("Sunz correction not done by satpy reader for versions >= 0.60.")

N_DETECTORS = 24

GEOLOCATION_NAMES = ['satellite_zenith_angle',
                     'solar_zenith_angle',
                     'satellite_azimuth_angle',
                     'solar_azimuth_angle',
                     'lat_pixels',
                     'lon_pixels']

PPS_TAGS = {"vii_668": "ch_r06",
            "vii_865": "ch_r09",
            "vii_1375": "ch_r13",
            "vii_1630": "ch_r16",
            "vii_3740": "ch_tb37",
            "vii_8540": "ch_tb85",
            "vii_10690": "ch_tb11",
            "vii_12020": "ch_tb12",
            "vii_2250": "ch_r22",
            # Not used yet:
            "vii_6725": "ch_tb67",
            "vii_13345": "ch_tb133",
            "vii_7325": "ch_tb73",
            "vii_443": "ch_rxx",
            "vii_555": "ch_rxx",
            "vii_752": "ch_rxx",
            "vii_763": "ch_rx",
            "vii_914": "ch_rxx",
            "vii_1240": "ch_rxx",
            "vii_3959": "ch_tbxx",
            "vii_4050": "ch_tbxx"}
IR_BANDS = ["vii_3740", "vii_8540", "vii_10690", "vii_12020", "vii_6725",
            "vii_13345", "vii_7325", "vii_3959", "vii_4050"]
refl_bands = get_refl_bands(PPS_TAGS)
ONE_IR_CHANNEL = "vii_10690"


def destripe(scene, band, n_scans_per_block=2):
    """
    Destripe data from given IR band using a detector mean filter.
    Filtering over two times number of detectors seem to give best results,
    as every second scan has higher offsets.
    """
    logger.info(f"Destriping IR channel {band}.")
    arr = scene[band].values
    n_rows = N_DETECTORS * n_scans_per_block
    n_blocks = int(np.floor(arr.shape[0] / n_rows))
    mean_row = np.nanmean(arr, axis=1)
    mean_per_row_in_block = np.nanmean(mean_row[0:n_blocks * n_rows].reshape(n_blocks, n_rows), axis=0)
    offset_per_row_in_block = mean_per_row_in_block - np.nanmean(mean_per_row_in_block)
    offset_per_pixel = np.tile(offset_per_row_in_block, (arr.shape[1], n_blocks + 1)).T[:arr.shape[0]]
    scene[band].values -= offset_per_pixel
    return offset_per_row_in_block


def set_header_and_band_attrs(scene, orbit_n=00000):
    """Set and delete some attributes."""
    nimg = 0  # name of first dataset is image0
    # Set some header attributes:
    irch = scene[ONE_IR_CHANNEL]
    scene.attrs['source'] = "metimage2pps.py"
    nimg = set_header_and_band_attrs_defaults(scene, PPS_TAGS, irch, orbit_n=orbit_n)
    for band in refl_bands:
        if band not in scene:
            continue
        # Sunz correction not done by satpy reader for versions >= 0.60.
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
    return nimg


def load_data(scene_files, channel_selection):
    """Load data."""
    scene = Scene(reader='vii_l1b_nc', filenames=scene_files)
    my_bands = get_band_names(PPS_TAGS, channel_selection)
    bands_to_load = my_bands + GEOLOCATION_NAMES
    scene.load(bands_to_load)
    return scene


def process_one_scene(scene_files, out_path,
                      engine='h5netcdf',
                      channel_selection="default",
                      destripe_ir_channels=False,
                      orbit_n=0,
                      platform_name=None):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    check_file_exists(scene_files)
    scene = load_data(scene_files, channel_selection)
    # one ir channel
    irch = scene[ONE_IR_CHANNEL]
    set_header_and_band_attrs(scene, orbit_n=orbit_n)
    rename_latitude_longitude(scene)
    adjust_lons_to_valid_range(scene)
    convert_angles(scene, delete_azimuth=True)
    update_angle_attributes(scene, irch)
    if destripe_ir_channels:
        for channel in IR_BANDS:
            if channel in scene:
                destripe(scene, channel)
    apply_sunz_correction(scene, refl_bands)
    if platform_name is not None:
        scene.attrs['platform'] = platform_name
    filename = compose_filename(scene, out_path, instrument='metimage', band=irch)
    header_attrs = get_header_attrs(scene, band=irch, sensor='metimage')
    save_data(scene, filename, header_attrs, engine)
    log_time(filename, tic)
    return filename
