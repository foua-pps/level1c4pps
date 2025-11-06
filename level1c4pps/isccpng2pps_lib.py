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
from pyorbital.astronomy import get_alt_az, sun_zenith_angle

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

channel_name = {"refl_00_65um": "VIS006",
                "refl_00_86um": "VIS008",
                "refl_01_60um": "IR_016"}

platform_id = {55: 321,
               70: 324}
# Normally read from file (HRIT). We do not have access to nominal calibration
# So we need to have that here. So far only 2021 
calibration_nominal = {2021: {321: {"VIS006": 24.974,
                                    "VIS008": 32.377,
                                    "IR_016": 23.710},
                              324: {"VIS006": 21.241,
                                    "VIS008": 27.921,
                                    "IR_016": 23.112}
                              }}
    
coef_slope_chan = ['refl_00_65um', 'refl_00_86um', 'refl_01_60um',
                   'temp_03_80um', 'temp_06_20um', 'temp_07_30um', 'temp_08_60um',
                   'temp_09_70um', 'temp_11_00um', 'temp_12_00um', 'temp_13_30um']
coef_slope_list = {
    '270_day': [1.00152, 0.946758, 1.00669, 0.926764, 0.991329, 0.999300, 1.01977, 1.03453, 0.996781, 1.02343, 0.863748],
    '270_nig': [1.00000, 1.00000, 1.00000, 0.950274, 0.989769, 1.00713, 1.01844, 1.01421, 0.995356, 1.02101, 0.905866],
    '271_day': [1.00298, 0.946479, 1.00762, 0.928358, 0.991436, 0.999576, 1.01939, 1.03306, 0.996656, 1.02405, 0.894260],
    '271_nig': [1.00000, 1.00000, 1.00000, 0.950015, 0.989834, 1.00698, 1.01804, 1.01361, 0.995194, 1.02155, 0.930947],
    '173_day': [1.00164, 0.944191, 1.00680, 0.929371, 1.00905, 1.00090, 1.01065, 1.02516, 0.996015, 1.02885, 0.871186],
    '173_nig': [1.00000, 1.00000, 1.00000, 0.950006, 1.00397, 1.00545, 1.00990, 1.01022, 0.994367, 1.02563, 0.912099],
    '55_day': [0.999586, 0.999657, 0.997100, 1.01192, 1.00193, 1.00093, 1.00202, 0.998663, 0.999412, 0.999663, 0.945356],
    '55_nig': [1.00000, 1.00000, 1.00000, 1.00772, 1.00129, 1.00127, 1.00180, 0.999506, 0.999569, 0.999708, 0.961137],
}
coef_offset_list = {
    '270_day': [0.000369844, 0.193868, 0.708022, 18.4795, 1.58729, -0.0312020, -4.93224, -7.83331, 1.13395, -5.51940, 30.5574],
    '270_nig': [0.00000, 0.00000, 0.00000, 11.8884, 1.97083, -2.13396, -4.75410, -2.48929, 1.46463, -4.95287, 19.1043],
    '271_day': [-0.00818343, 0.194117, 0.718540, 17.7929, 1.57229, -0.0877037, -4.83326, -7.49875, 1.17318, -5.66102, 23.6271],
    '271_nig': [0.00000, 0.00000, 0.00000, 11.9196, 1.96629, -2.07708, -4.64315, -2.38190, 1.51330, -5.07539, 13.6827],
    '173_day': [-4.45259e-05, 0.214025, 0.703233, 17.3952, -2.21718, -0.347357, -2.62520, -5.68860, 1.37321, -6.76384, 28.8614],
    '173_nig': [0.00000, 0.00000, 0.00000, 11.9012, -0.963594, -1.56746, -2.52102, -1.75993, 1.76101, -5.99385, 17.7467],
    '55_day': [0.00410180, 0.0160834, 0.101325, -2.83943, -0.417072, -0.181784, -0.486767, 0.304194, 0.134653, 0.0706959, 12.3763],
    '55_nig': [0.00000, 0.00000, 0.00000, -1.61703, -0.258592, -0.278584, -0.442302, 0.0806762, 0.0877965, 0.0597462, 8.17879]
}

coef_slope_list_v1 = {
    '270_day': [1.00166, 0.948021, 1.00000, 0.926764, 0.991329, 0.999300,
                1.01977, 1.03453, 0.996781, 1.02343, 0.863748],
    '270_all': [1.00000, 1.00000, 1.00000, 0.940848, 0.990431, 1.00483,
                1.01947, 1.01835, 0.996272, 1.02239, 0.890045],
    '271_day': [1.00313, 0.947739, 1.00000, 0.928358, 0.991436, 0.999576,
                1.01939, 1.03306, 0.996656, 1.02405, 0.894260],
    '271_all': [1.00000, 1.00000, 1.00000, 0.939285, 0.990504, 1.00481,
                1.01907, 1.01759, 0.996129, 1.02298, 0.917607],
    '173_day': [1.00169, 0.945509, 1.00000, 0.929371, 1.00905, 1.00090,
                1.01065, 1.02516, 0.996015, 1.02885, 0.871186],
    '173_all': [1.00000, 1.00000, 1.00000, 0.938706, 1.00526, 1.00413, 1.01050,
                1.01338, 0.995399, 1.02742, 0.896835],
    '55_day': [0.999577, 0.999643, 1.00000, 1.01192, 1.00193, 1.00093, 1.00202,
               0.998663, 0.999412, 0.999663, 0.945356],
    '55_all': [1.00000, 1.00000, 1.00000, 1.00882, 1.00146, 1.00121, 1.00194,
               0.999321, 0.999513, 0.999684, 0.955400]
}
coef_offset_list_v1 = {
    '270_day': [-0.00207069, 0.161403, 0.00000, 18.4795, 1.58729, -0.0312020,
                -4.93224, -7.83331, 1.13395, -5.51940, 30.5574],
    '270_all': [0.00000, 0.00000, 0.00000, 14.4615, 1.80946, -1.50406,
                -4.94939, -3.61878, 1.23287, -5.28908, 23.4361],
    '271_day': [-0.0107727, 0.161054, 0.00000, 17.7929, 1.57229, -0.0877037,
                -4.83326, -7.49875, 1.17318, -5.66102, 23.6271],
    '271_all': [0.00000, 0.00000, 0.00000, 14.7561, 1.80263, -1.48189,
                -4.84116, -3.46801, 1.27653, -5.42165, 17.3323],
    '173_day': [-0.00100716, 0.162379, 0.00000, 17.3952, -2.21718, -0.347357,
                -2.62520, -5.68860, 1.37321, -6.76384, 28.8614],
    '173_all': [0.00000, 0.00000, 0.00000, 14.8503, -1.29354, -1.20664,
                -2.63734, -2.61762, 1.49858, -6.43452, 21.9239],
    '55_day': [0.00444712, 0.0162967, 0.00000, -2.83943, -0.417072, -0.181784,
               -0.486767, 0.304194, 0.134653, 0.0706959, 12.3763],
    '55_all': [0.00000, 0.00000, 0.00000, -1.94343, -0.302762, -0.259928,
               -0.471424, 0.131567, 0.104622, 0.0658078, 9.71430]
}
satellite_names = {270: "GOES-16",  # ABI
                   271: "GOES-17",  # ABI
                   173: "Himawari-8",  # AHI
                   55: "Meteosat-8",  # MSG1 SEVIRI
                   70: "Meteosat-11"}


def get_encoding_isccpng(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene, BANDNAMES, PPS_TAGNAMES, chunks=None)


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
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
    return nimg


def homogenize_channel(scene, wmo_id, illum, band, sol_zen):
    """Homogenize one satellite for one band."""
    data = scene[band]
    dwmo_id = scene["wmo_id"].values
    ch_index = coef_slope_chan.index(band)
    k = coef_slope_list[str(wmo_id) + illum][ch_index]
    c = coef_offset_list[str(wmo_id) + illum][ch_index]
    logger.info(
        f"Homogenizing channel {band} for {satellite_names[wmo_id]} for illum:{illum}"
        f" with respect to SEVIRI on MSG4 using y={k} * x + {c}")
    dn_discr = 90
    if illum == '_day':
        update = (dwmo_id == wmo_id) & (data > 0) & (sol_zen < dn_discr)
        scene[band].values = np.where(update, scene[band].values * k + c, scene[band].values)

    if illum == '_nig':
        update = (dwmo_id == wmo_id) & (data > 0) & (sol_zen >= dn_discr)
        scene[band].values = np.where(update, scene[band].values * k + c, scene[band].values)


def homogenize(scene):
    """Homogenize data to Meteosat-11."""
    sol_zen = scene["sunzenith"]
    for band in BANDNAMES:
        if band not in scene:
            continue
        for wmo_id in satellite_names:
            if wmo_id == 70:
                continue  # Meteosat11
            for illumination in ["_day", "_nig"]:
                homogenize_channel(scene, wmo_id, illumination, band, sol_zen)


def recalibrate_meteosat(scene):
    """Nominal calibration is applied, redo with meirnik calibration."""
    # from satpy.readers.seviri_base import MeirinkCoefficients
    from satpy.readers.core.seviri import MeirinkCoefficients
    start_time = dt64_to_datetime(scene["refl_00_65um"].attrs["start_time"])
    for wmo_id in platform_id:
        for band in channel_name:
            if band not in scene:
                continue
            dwmo_id = scene["wmo_id"]
            old_gain = calibration_nominal[start_time.year][platform_id[wmo_id]][channel_name[band]]
            meirink = MeirinkCoefficients(platform_id[wmo_id], channel_name[band], start_time)
            new_gain = 1000.0 * meirink._get_gain()['MEIRINK-2023']
            update = scene["wmo_id"].values == wmo_id
            logger.info(
                f"Recalibrating channel {band} for {satellite_names[wmo_id]} using y *= {new_gain} /{old_gain}")
            scene[band].values = np.where(update, scene[band].values * new_gain / old_gain, scene[band].values)


def get_solar_angles(scene, lons, lats):
    """Compute solar angles.
    Compute angles for each scanline using their acquisition time to account for
    the earth's rotation over the course of one scan.
    Returns:
        Solar azimuth angle, Solar zenith angle in degrees
    """
    acq_time =  scene["pixel_time"].copy()
    _, suna = get_alt_az(acq_time, lons, lats)
    suna = np.rad2deg(suna)
    sunz = sun_zenith_angle(acq_time, lons, lats)
    return suna, sunz



def fix_pixel_time(scene):
    """Fix the time pixel variable, original file does not contain units."""
    del scene["pixel_time"].coords["crs"]
    scene["pixel_time"].encoding['coordinates'] = "lon lat"
    scene["pixel_time"].data = scene["pixel_time"].data * np.timedelta64(1, 's') + scene['temp_11_00um'].attrs["start_time"]

    
def update_solar_angles(scene):
    """USe pixel time to calculate solar angles."""
    suna, sunz = get_solar_angles(scene, lons=scene["lon"], lats=scene["lat"])
    scene["solar_zenith_angle"].values = sunz.values
    scene["solar_azimuth_angle"].values = suna.values


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
    update_solar_angles(scn_)
    convert_angles(scn_, delete_azimuth=True)
    update_angle_attributes(scn_, irch)
    recalibrate_meteosat(scn_)
    homogenize(scn_)
    apply_sunz_correction(scn_, REFL_BANDS)

    filename = compose_filename(scn_, out_path, instrument='seviri', band=irch)
    encoding = get_encoding_isccpng(scn_)

    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='seviri'),
                       engine=engine,
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=encoding)
    logging.info("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
    return filename
