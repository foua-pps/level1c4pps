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

"""Functions to convert VIIRS level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import time
from satpy.scene import Scene
from level1c4pps import (compose_filename,
                         set_header_and_band_attrs_defaults,
                         rename_latitude_longitude,
                         save_data,
                         get_refl_bands,
                         log_time,
                         get_band_names,
                         update_angle_attributes, get_header_attrs,
                         convert_angles)
import pyspectral  # testing that pyspectral is available # noqa: F401
import logging

# Example:

logger = logging.getLogger('viirs2pps')


ANGLE_NAMES = ['satellite_zenith_angle',
               'solar_zenith_angle',
               'satellite_azimuth_angle',
               'solar_azimuth_angle']

PPS_TAGNAMES = {"M05": 'ch_r06',
                "M07": 'ch_r09',
                "M09": 'ch_r13',
                "M10": 'ch_r16',
                "M12": 'ch_tb37',
                "M11": 'ch_r22',
                "M14": 'ch_tb85',
                "M15": 'ch_tb11',
                "M16": 'ch_tb12',
                "I01": 'ch_r06',
                "I02": 'ch_r09',
                "I03": 'ch_r16',
                "I04": 'ch_tb37',
                # Not used by pps:
                "I05": 'ch_tbxx',
                "M01": 'ch_rxx',
                "M02": 'ch_rxx',
                "M03": 'ch_rxx',
                "M04": 'ch_rxx',
                "M06": 'ch_rxx',
                "M08": 'ch_rxx',
                "M13": 'ch_tbxx'}
refl_bands = get_refl_bands(PPS_TAGNAMES)
ONE_IR_CHANNEL = 'M15'


def set_header_and_band_attrs(scene, orbit_n=0):
    """Set and delete some attributes."""
    irch = scene[ONE_IR_CHANNEL]
    nimg = set_header_and_band_attrs_defaults(scene, PPS_TAGNAMES, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "viirs2pps.py"
    if 'I04' in scene:
        # If highresolution we should have I04,
        scene.attrs['number_of_scans'] = scene['I04'].values.shape[0] / scene['I04'].attrs['rows_per_scan']
    else:
        # else use 11um.
        scene.attrs['number_of_scans'] = scene["M15"].values.shape[0] / scene["M15"].attrs['rows_per_scan']
    for band in refl_bands:
        if band not in scene:
            continue
        # VIIRS data read with sunz_corrected modifier
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'True'
    return nimg


def load_data(scene_files, reader, all_channels, pps_channels, use_iband_res):
    """Load data."""
    scene = Scene(
        reader=reader,
        filenames=scene_files)
    my_bands = get_band_names(PPS_TAGNAMES, all_channels, pps_channels)
    if use_iband_res:
        my_ibands_i = [band for band in my_bands if "I" in band]
        my_ibands_i_tags = [PPS_TAGNAMES[band] for band in my_ibands_i]
        my_ibands_m = [band for band in my_bands if "M" in band and PPS_TAGNAMES[band] not in my_ibands_i_tags]
    else:
        my_bands = [band for band in my_bands if "I" not in band]

    if use_iband_res:
        my_ibands_to_load = my_ibands_i + ANGLE_NAMES + ['i_latitude', 'i_longitude']
        scene.load(my_ibands_to_load, resolution=371)
        scene.load(my_ibands_m, resolution=742)
        scene = scene.resample(resampler='native')
    elif reader == "viirs_compact":
        my_mband_refl = [band for band in my_bands if band in refl_bands]
        my_mband_tb = [band for band in my_bands if band not in refl_bands]
        my_unmodified_bands_to_load = my_mband_tb + ANGLE_NAMES + ['latitude_m', 'longitude_m']
        scene.load(my_unmodified_bands_to_load, resolution=742)
        # Load reflective bands with sunz-correction (not the default for VIIRS compact).
        scene.load(my_mband_refl, modifiers=("sunz_corrected", ))
    else:
        bands_to_load = my_bands + ANGLE_NAMES + ['m_latitude', 'm_longitude']
        scene.load(bands_to_load, resolution=742)
    return scene


def process_one_scene(scene_files, out_path, use_iband_res=False, reader='viirs_sdr', engine='h5netcdf',
                      all_channels=False, pps_channels=False, orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scene = load_data(scene_files, reader, all_channels, pps_channels, use_iband_res)
    irch = scene[ONE_IR_CHANNEL]
    set_header_and_band_attrs(scene, orbit_n=orbit_n)
    rename_latitude_longitude(scene)
    convert_angles(scene, delete_azimuth=True)
    update_angle_attributes(scene, irch)
    header_attrs = get_header_attrs(scene, band=irch, sensor='viirs')
    filename = compose_filename(scene, out_path, instrument='viirs', band=irch)
    save_data(scene, filename, header_attrs, engine)
    log_time(filename, tic)
    return filename
