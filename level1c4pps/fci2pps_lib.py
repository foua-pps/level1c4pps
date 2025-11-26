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

"""Functions to convert FCI level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import os
import time
import satpy
from satpy.scene import Scene
from level1c4pps import (get_encoding,
                         compose_filename,
                         get_header_attrs,
                         set_header_and_band_attrs_defaults,
                         dt64_to_datetime,
                         fix_timestamp_datatype)
from level1c4pps.seviri2pps_lib import (get_lonlats,
                                        add_ancillary_datasets,
                                        get_satellite_angles,
                                        make_azidiff_angle)
import logging
import numpy as np
from pyorbital.astronomy import get_alt_az, sun_zenith_angle
import hdf5plugin  # testing that library for fci is available # noqa: F401
import datetime as dt
import pytz  # testing that library for fci is available # noqa: F401
from packaging.version import Version

# from satpy.utils import debug_on
# debug_on()

# Example:
# 'W_XX-EUMETSAT-Darmstadt,IMG+SAT,MTI1+FCI-1C-RRAD-HRFI-FD--CHK-BODY--DIS-NC4E_C_EUMT_20250408093037_IDPFI_OPE_20250408092752_20250408092842_N_JLS_O_0057_0034.nc"


logger = logging.getLogger('fci2pps')

if Version(satpy.__version__) < Version('0.59.1'):
    if Version(satpy.__version__) > Version('0.56.0'):
        logger.warning("Native resampling craching for satpy 0.57 to 0.59.0.")


BANDNAMES_DEFAULT = ["vis_06",
                     "vis_08",
                     "nir_13",
                     "nir_16",
                     "ir_38",
                     "wv_73",

                     "ir_105",
                     "ir_123",
                     "nir_22",
                     "wv_63",
                     "ir_133"]

BANDNAMES_PPS = ["vis_06",
                 "vis_08",
                 "nir_13",
                 "nir_16",
                 #"nir_22",
                 "ir_38",
                 "ir_87",
                 "ir_105",
                 "ir_123"]

REFL_BANDS = ["vis_06", "vis_08", "nir_13", "nir_16", "vis_443",
              "vis_05", "vis_04", "vis_09",
              "nir_22"]

ANGLE_NAMES = ['observation_zenith', 'solar_zenith',
               'observation_azimuth', 'solar_azimuth']

PPS_TAGNAMES = {"vis_06": "ch_r06",
                "vis_08": "ch_r09",
                "nir_13": "ch_r13",
                "nir_16": "ch_r16",
                "ir_38": "ch_tb37",
                "ir_87": "ch_tb85",
                "ir_105": "ch_tb11",
                "ir_123": "ch_tb12",
                "nir_22": "ch_r22",
                # Not used yet:
                "wv_63": "ch_tb67",
                "ir_133": "ch_tb133",
                "wv_73": "ch_tb73",
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
    nimg = 0  # name of first dataset is image0
    # Set some header attributes:
    irch = scene['ir_105']
    scene.attrs['source'] = "fci2pps.py"
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    for band in REFL_BANDS:
        if band not in scene:
            continue
        print("Is this correct???????????????")
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'True'
    return nimg


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
    if type(scene["ir_105_time"].values[0, 0]) == np.float64:
        epoch_to_2000 = dt.datetime(2000, 1, 1, tzinfo=dt.timezone.utc).timestamp()
        scene["ir_105_time"] = scene["ir_105_time"] + epoch_to_2000
        scene["ir_105_time"] = scene["ir_105_time"].astype('datetime64[s]')
    ind = np.int16(scene["ir_105_time"].shape[0]/2)
    a_time = dt64_to_datetime(scene["ir_105_time"].values[ind, ind])
    if a_time > scene.end_time or a_time < scene.start_time:
        raise ValueError
    scene["ir_105_time"].attrs.pop("units", None)


def resample_data(scn_in, datasets, resample_grid="coarse"):
    logger.info("Resampling to coarsest area")
    if resample_grid in ["fine"]:
        logger.info("Resampling to finest grid")
        scn_out = scn_in.resample(scn_in.finest_area(), datsets=datasets, resampler='native')
    elif resample_grid in ["coarse"]:
        logger.info("Resampling to coarsest grid")
        scn_out = scn_in.resample(scn_in.coarsest_area(), datsets=datasets, resampler='native')
    elif "msg" in resample_grid:
        import psutil
        if psutil.virtual_memory()[0]/32e9 < 1:
            logger.warning(
                "Resampling to msg grid, via native resampling to coarsest grid to save RAM. "
                "As we have less than 32G RAM.")
            scn_out = scn_in.resample(scn_in.coarsest_area(), datsets=datasets, resampler='native')
            scn_out = scn_out.resample(resample_grid, datsets=datasets, resampler='nearest')
        else:
            logger.info("Resampling to msg grid")
            scn_out = scn_in.resample(resample_grid, datsets=datasets, resampler='nearest')
    return scn_out


def process_one_scene(scene_files, out_path,
                      engine='h5netcdf',
                      all_channels=False,
                      pps_channels=False,
                      resample_grid=False,
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
    scn_ = resample_data(scn_in, MY_BANDNAMES + ["ir_105_time"], resample_grid=resample_grid)
    fix_time(scn_)
    irch = scn_['ir_105']
    lons, lats = get_lonlats(scn_['ir_105'])
    suna, sunz = get_solar_angles(scn_, lons=lons, lats=lats)
    sata, satz = get_satellite_angles(scn_['ir_105'], lons=lons, lats=lats)
    azidiff = make_azidiff_angle(sata, suna)
    sata = None
    suna = None
    add_ancillary_datasets(scn_,
                           lons=lons, lats=lats,
                           sunz=sunz, satz=satz,
                           azidiff=azidiff,
                           suna=suna, sata=sata,
                           irch_name="ir_105",
                           save_azimuth_angles=False,
                           chunks=(464, 928))
    set_header_and_band_attrs(scn_, orbit_n=orbit_n)
    scn_["ir_105_time"].attrs.pop("_FillValue", None)
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
