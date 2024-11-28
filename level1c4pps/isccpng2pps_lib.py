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

#   Nina Hakansson <nina.hakansson@smhi.se>

"""Functions to convert ISCCP Next Generation level-1-G data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import os
import time
import xarray as xr
import numpy as np
from satpy.scene import Scene
from level1c4pps import (get_encoding,
                         dt64_to_datetime,
                         compose_filename,
                         apply_sunz_correction,
                         rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         set_header_and_band_attrs_defaults,
                         convert_angles,
                         adjust_lons_to_valid_range)

import logging

# from satpy.utils import debug_on
# debug_on()

# Example:
# 'SEVIRI-like_ISCCP-NG_L1g_demo_v2_res_0_05deg__temp_13_30um__20200701T1200.nc'


logger = logging.getLogger('isccpng2pps')

BANDNAMES = ['refl_01_60um',
             'refl_00_65um',
             'refl_00_86um',
             'temp_03_80um',
             'temp_08_60um',
             'temp_09_70um',
             'temp_11_00um',
             'temp_12_00um',
             'temp_13_30um',                    
             'temp_06_20um',
             'temp_07_30um']


REFL_BANDS = ['refl_01_60um',
              'refl_00_65um',
              'refl_00_86um']

ANGLE_NAMES = ['solar_zenith_angle', 'satellite_zenith_angle',
               'solar_azimuth_angle', 'satellite_azimuth_angle']

PPS_TAGNAMES = {'refl_01_60um': "ch_r16",
                'refl_00_65um': "ch_r06",
                'refl_00_86um': "ch_r09",
                'temp_03_80um': "ch_tb37",
                'temp_08_60um': "ch_tb85",
                'temp_09_70um': "ch_tbxx",
                'temp_11_00um': "ch_tb11",
                'temp_12_00um': "ch_tb12",
                'temp_13_30um': "ch_tb133",                    
                'temp_06_20um': "ch_tbxx",
                'temp_07_30um': "ch_tb73"}

BANDNAMES = list(PPS_TAGNAMES.keys())


def get_encoding_isccpng(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene, orbit_n=00000):
    """Set and delete some attributes."""
    nimg = 0  # name of first dataset is image0
    # Set some header attributes:
    irch = scene['temp_11_00um']
    irch.attrs['instrument'] = "seviri"
    scene.attrs['source'] = "isccpng2pps.py"
    scene.attrs['platform_name'] = "meteosat11"
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    for band in REFL_BANDS:
        if band not in scene:
            continue
        print("Is this correct, we think so")
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
    return nimg

def homogenize(scene):
    pass

def recalibrate_meteosat(scene):
    """nominal calibration is applied, redo with meirnik calibration."""
    from satpy.readers.seviri_base import MeirinkCalibrationHandler
    start_time = dt64_to_datetime(scene["refl_00_65um"].attrs["start_time"])
    channel_name = {"refl_00_65um": "VIS006",
                    "refl_00_86um": "VIS008",
                    "refl_01_60um": "IR_016"}
    platform_id = {55: 321,
                   70: 324}
    # Normally read from file (HRIT). We do not have access to nominal calibration
    calibration_nominal = {2021: {321: {"VIS006": 24.974,
                                        "VIS008": 32.377,
                                        "IR_016": 23.710},
                                  324: {"VIS006": 21.241,  
                                        "VIS008": 27.921,
                                        "IR_016": 23.112}
                                   }}
    
    for wmo_id in platform_id:        
        for band in channel_name:
            old_gain = calibration_nominal[start_time.year][platform_id[wmo_id]][channel_name[band]]
            meirink = MeirinkCalibrationHandler(calib_mode="MEIRINK-2023")
            new_gain = 1000.0 * meirink.get_slope(platform_id[wmo_id], channel_name[band], start_time)
            scene[band] = scene[band] * new_gain / old_gain



def fix_pixel_time(scene):
    """Fix the time pixel variable, original file does not contain units."""
    del scene["pixel_time"].coords["crs"]
    scene["pixel_time"].encoding['coordinates'] = "lon lat"
    scene["pixel_time"].data = scene["pixel_time"].data * np.timedelta64(1, 's') + scene['temp_11_00um'].attrs["start_time"]
    
def process_one_scene(scene_files, out_path,
                      engine='h5netcdf',
                      orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(reader='multiple_sensors_isccpng_l1g_nc', filenames=scene_files)
    scn_.load(BANDNAMES + ANGLE_NAMES + ["wmo_id", "pixel_time", "lon", "lat"])

    # one ir channel
    irch = scn_['temp_11_00um']

    set_header_and_band_attrs(scn_, orbit_n=orbit_n)
    fix_pixel_time(scn_)
    # rename_latitude_longitude(scn_)
    adjust_lons_to_valid_range(scn_)
    convert_angles(scn_, delete_azimuth=True)
    update_angle_attributes(scn_, irch)
    apply_sunz_correction(scn_, REFL_BANDS)
    recalibrate_meteosat(scn_)
    homogenize(scn_)
    
    filename = compose_filename(scn_, out_path, instrument='seviri', band=irch)
    encoding = get_encoding_isccpng(scn_)

    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='seviri'),
                       engine=engine,
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=encoding)
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
    return filename
