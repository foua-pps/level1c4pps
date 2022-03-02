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
                         rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         convert_angles)
import pyspectral  # testing that pyspectral is available # noqa: F401
import logging

# Example:

logger = logging.getLogger('viirs2pps')

# Order of BANDNAMES decides order of channels in file. Not important
# but nice to have the same order for I- and M-bands
BANDNAMES = ["M01",  "M02", "M03", "M04",
             "M05", "M06", "M07",  # 0.6, 0.7, 0.9 M-band
             "I01", "I02",         # 0.6, 0.9 I-band
             "M08", "M09",         # 1.2, 1.3 M-band
             "M10",                # 1.6 M-band
             "I03",                # 1.6 I-band
             "M11",                # 2.25 M-band
             "M12",                # 3.7 M-band
             "I04",                # 3.7 I-band
             "M13", "M14",         # 4.05, 8.55 M-band
             "M15", "M16",         # 11, 12 M-band
             "I05"]                # 11.5 I-band

MBANDS = ["M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08",
          "M09", "M10", "M11", "M12", "M13", "M14", "M15", "M16"]

IBANDS = ["I01", "I02", "I03", "I04", "I05"]

REFL_BANDS = ["M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08",
              "M09", "M10", "M11", "I01", "I02", "I03"]

MBAND_PPS = ["M05", "M07", "M09", "M10", "M12", "M14", "M15", "M16"]
IBAND_PPS_I = ["I01", "I02", "I03", "I04"]
IBAND_PPS_M = ["M09", "M14", "M15", "M16"]

MBAND_DEFAULT = ["M05", "M07", "M09", "M10",  "M11", "M12", "M14", "M15", "M16"]
IBAND_DEFAULT_I = ["I01", "I02", "I03", "I04"]
IBAND_DEFAULT_M = ["M09",   "M11", "M14", "M15", "M16"]

ANGLE_NAMES = ['satellite_zenith_angle', 'solar_zenith_angle',
               'satellite_azimuth_angle', 'solar_azimuth_angle']

PPS_TAGNAMES = {"M05": 'ch_r06',
                "M07": 'ch_r09',
                "M09": 'ch_r13',
                "M10": 'ch_r16',
                "M12": 'ch_tb37',
                "M14": 'ch_tb85',
                "M15": 'ch_tb11',
                "M16": 'ch_tb12',
                "I01": 'ch_r06',
                "I02": 'ch_r09',
                "I03": 'ch_r16',
                "I04": 'ch_tb37',
                # Not used by pps:
                "M11": 'ch_r22',
                "I05": 'ch_tbxx',
                "M01": 'ch_rxx',
                "M02": 'ch_rxx',
                "M03": 'ch_rxx',
                "M04": 'ch_rxx',
                "M06": 'ch_rxx',
                "M08": 'ch_rxx',
                "M13": 'ch_tbxx'}


def get_encoding_viirs(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene, orbit_n=0):
    """Set and delete some attributes."""
    irch = scene['M15']
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "viirs2pps.py"
    if 'I04' in scene:
        # If highresolution we should have I04,
        scene.attrs['number_of_scans'] = scene['I04'].values.shape[0]/scene['I04'].attrs['rows_per_scan']
    else:
        # else use 11um.
        scene.attrs['number_of_scans'] = scene['M15'].values.shape[0]/scene['M15'].attrs['rows_per_scan']
    for band in REFL_BANDS:
        if band not in scene:
            continue
        # For VIIRS data sun_zenith_angle_correction_applied is applied always!
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'True'
    return nimg


def process_one_scene(scene_files, out_path, use_iband_res=False, engine='h5netcdf',
                      all_channels=False, pps_channels=False, orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(
        reader='viirs_sdr',
        filenames=scene_files)

    MY_MBAND = MBAND_DEFAULT
    MY_IBAND_I = IBAND_DEFAULT_I
    MY_IBAND_M = IBAND_DEFAULT_M

    if all_channels:
        MY_MBAND = MBANDS
        MY_IBAND_I = IBANDS
        MY_IBAND_M = MBANDS
    if pps_channels:
        MY_MBAND = MBAND_PPS
        MY_IBAND_I = IBAND_PPS_I
        MY_IBAND_M = IBAND_PPS_M

    if use_iband_res:
        scn_.load(MY_IBAND_I + ANGLE_NAMES + ['i_latitude', 'i_longitude'], resolution=371)
        scn_.load(MY_IBAND_M, resolution=742)
        scn_ = scn_.resample(resampler='native')
    else:
        scn_.load(MY_MBAND + ANGLE_NAMES + ['m_latitude', 'm_longitude'], resolution=742)

    # one ir channel
    irch = scn_['M15']

    # Set header and band attributes
    set_header_and_band_attrs(scn_, orbit_n=orbit_n)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, delete_azimuth=True)
    update_angle_attributes(scn_, irch)

    filename = compose_filename(scn_, out_path, instrument='viirs', band=irch)
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='viirs'),
                       engine=engine,
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_viirs(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
    return filename
