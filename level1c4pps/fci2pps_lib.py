#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2025 level1c4pps developers
#
# This file is part of level1c4pps
#
# level1c4pps is free software: you can redistribute it and/or modify it
# under the terms of the Apache-2.0 license.
#
#

"""Functions to convert FCI level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import datetime as dt
import logging
import os
import time

import hdf5plugin  # testing that library for fci is available # noqa: F401
import numpy as np
import pytz  # testing that library for fci is available # noqa: F401
import satpy
from packaging.version import Version
from pyorbital.astronomy import get_alt_az, sun_zenith_angle
from satpy.scene import Scene

from level1c4pps import (compose_filename, dt64_to_datetime,
                         fix_timestamp_datatype, get_encoding,
                         get_header_attrs, set_header_and_band_attrs_defaults)
from level1c4pps.seviri2pps_lib import (add_ancillary_datasets, get_lonlats,
                                        get_satellite_angles,
                                        make_azidiff_angle)

# Example:
# 'W_XX-EUMETSAT-Darmstadt,IMG+SAT,MTI1+FCI-1C-RRAD-HRFI-FD--CHK-BODY--DIS-NC4E_C_EUMT_20250408093037_IDPFI_OPE_20250408092752_20250408092842_N_JLS_O_0057_0034.nc"

logger = logging.getLogger('fci2pps')

if Version(satpy.__version__) <= Version('0.59.0'):
    if Version(satpy.__version__) > Version('0.56.0'):
        logger.warning("Native resampling craching for satpy 0.57 to 0.59.0.")


BANDNAMES_DEFAULT = ["vis_06",
                     "vis_08",
                     "nir_13",
                     "nir_16",
                     "nir_22",
                     "ir_38",
                     "wv_63",
                     "wv_73",
                     "ir_105",
                     "ir_123",
                     "ir_133"]

BANDNAMES_PPS = ["vis_06",
                 "vis_08",
                 "nir_13",
                 "nir_16",
                 "nir_22",
                 "ir_38",
                 "ir_87",
                 "ir_105",
                 "ir_123"]

REFL_BANDS = ["vis_06", "vis_08", "nir_13", "nir_16",
              "vis_05", "vis_04", "vis_09",
              "nir_22"]

ANGLE_NAMES = ['observation_zenith', 'solar_zenith',
               'observation_azimuth', 'solar_azimuth']

PPS_TAGNAMES = {"vis_06": "ch_r06",
                "vis_08": "ch_r09",
                "nir_13": "ch_r13",
                "nir_16": "ch_r16",
                "nir_22": "ch_r22",
                "ir_38": "ch_tb37",
                "ir_87": "ch_tb85",
                "ir_105": "ch_tb11",
                "ir_123": "ch_tb12",
                # Not used yet:
                "wv_63": "ch_tb67",
                "wv_73": "ch_tb73",
                "ir_133": "ch_tb133",
                "vis_05": "ch_rxx",
                "vis_04": "ch_rxx",
                "vis_09": "ch_rxx"}

BANDNAMES = list(PPS_TAGNAMES.keys())


def get_encoding_fci(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene, orbit_n=00000):
    """Set and delete some attributes."""
    # Set some header attributes:
    irch = scene['ir_105']
    scene.attrs['source'] = "fci2pps.py"
    set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    for band in REFL_BANDS:
        if band not in scene:
            continue
        if scene[band].attrs['sun_zenith_angle_correction_applied'] == "False":
            logger.warning(f"Setting sun_zenith_angle_correction_applied to True for band {band}")
            scene[band].attrs['sun_zenith_angle_correction_applied'] = 'True'


def get_solar_angles(scene, lons, lats):
    """Compute solar angles.

    Compute angles for each scanline using their acquisition time to account for
    the earth's rotation over the course of one scan.

    Returns:
        Solar azimuth angle, Solar zenith angle in degrees

    """
    suna = np.full(lons.shape, np.nan)
    sunz = np.full(lons.shape, np.nan)
    acq_time = scene["ir_105_time"].copy()
    _, suna = get_alt_az(acq_time, lons, lats)
    suna = np.rad2deg(suna)
    sunz = sun_zenith_angle(acq_time, lons, lats)
    return suna, sunz


def fix_time(scene):
    """Make datetime objects from time in seconds."""
    if type(scene["ir_105_time"].values[0, 0]) is np.float64:
        epoch_to_2000 = dt.datetime(2000, 1, 1, tzinfo=dt.timezone.utc).timestamp()
        scene["ir_105_time"] = (scene["ir_105_time"] + epoch_to_2000).astype('datetime64[s]')
    ind = int(scene["ir_105_time"].shape[0]/2)
    a_time = dt64_to_datetime(scene["ir_105_time"].values[ind, ind])
    if a_time > scene.end_time or a_time < scene.start_time:
        raise ValueError(f"There is someting wrong with the time variable, at {ind}, {ind}")
    scene["ir_105_time"].attrs.pop("units", None)
    scene["ir_105_time"].attrs.pop("_FillValue", None)


def resample_data(scn_in, datasets, resample_grid="coarse", resample_save_ram=False):
    """Resample data to the same grid."""
    if resample_grid not in ["coarse"]:
        logger.info(
            "More than 16GB of RAM memory is needed to resample to fine resolution.\n"
            "More than 16GB of RAM is needed to use nearset neighbour resampling directly to msg area.\n"
            "Try to run with --resample_via_native_coarse, to remap to msg area with less RAM usage.")

    if resample_grid in ["fine"]:
        logger.info("Resampling to finest grid")
        scn_out = scn_in.resample(scn_in.finest_area(), datsets=datasets, resampler="native")
    elif resample_grid in ["1km"]:
        logger.info("Resampling to 1km grid")
        scn_out = scn_in.resample(scn_in["vis_08"].area, datsets=datasets, resampler="native")
    elif resample_grid in ["coarse"]:
        logger.info("Resampling to coarsest grid")
        scn_out = scn_in.resample(scn_in.coarsest_area(), datsets=datasets, resampler="native")
    elif "msg" in resample_grid:
        if resample_save_ram:
            logger.info("Resampling to msg grid, via coarsest channel area.")
            scn_out = scn_in.resample(scn_in.coarsest_area(), datsets=datasets, resampler='native')
            scn_out = scn_out.resample(resample_grid, datsets=datasets, resampler='nearest')
        else:
            logger.info("Resampling to msg grid")
            scn_out = scn_in.resample(resample_grid, datsets=datasets, resampler='nearest')
    else:
        raise ValueError (f"No remapping method for grid {resample_grid}")
    return scn_out


def add_angles_and_latlon(scene):
    """Add the lon/lat and angles datasets."""
    lons, lats = get_lonlats(scene['ir_105'])
    suna, sunz = get_solar_angles(scene, lons=lons, lats=lats)
    sata, satz = get_satellite_angles(scene['ir_105'], lons=lons, lats=lats)
    azidiff = make_azidiff_angle(sata, suna)
    sata = None
    suna = None
    add_ancillary_datasets(scene,
                           lons=lons, lats=lats,
                           sunz=sunz, satz=satz,
                           azidiff=azidiff,
                           suna=suna, sata=sata,
                           irch_name="ir_105",
                           save_azimuth_angles=False,
                           chunks=(464, 928))


def process_one_scene(scene_files, out_path,
                      engine='h5netcdf',
                      all_channels=False,
                      pps_channels=False,
                      resample_grid="coarse",
                      resample_save_ram=False,
                      orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_in = Scene(reader='fci_l1c_nc', filenames=scene_files)
    MY_BANDNAMES = BANDNAMES_DEFAULT
    if all_channels:
        MY_BANDNAMES = BANDNAMES
    if pps_channels:
        MY_BANDNAMES = BANDNAMES_PPS
    scn_in.load(MY_BANDNAMES + ["ir_105_time"])
    scn_ = resample_data(scn_in, MY_BANDNAMES + ["ir_105_time"],
                         resample_grid=resample_grid,
                         resample_save_ram=resample_save_ram)
    fix_time(scn_)
    add_angles_and_latlon(scn_)
    irch = scn_['ir_105']
    set_header_and_band_attrs(scn_, orbit_n=orbit_n)
    filename = compose_filename(scn_, out_path, instrument='fci', band=irch)
    encoding = get_encoding_fci(scn_)
    fix_timestamp_datatype(scn_, encoding, "ir_105_time")
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='fci'),
                       engine=engine,
                       include_lonlats=False,
                       flatten_attrs=True,
                       pretty=True,
                       encoding=encoding)
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename), time.time() - tic))
    return filename
