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

"""Functions to convert VGAC level-1c data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

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

logger = logging.getLogger('vgac2pps')

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


REFL_BANDS = ["M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08",
              "M09"]

MBAND_PPS = ["M05", "M07", "M09", "M10", "M11", "M12", "M14", "M15", "M16"]

MBAND_AVHRR = ["M05", "M07", "M12", "M15", "M16"]

MBAND_DEFAULT = ["M05", "M07", "M09", "M10",  "M11", "M12", "M14", "M15", "M16"]


ANGLE_NAMES = ['vza', 'sza', 'azn', 'azi']

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

def convert_to_noaa19(scene):
    """
    Convert channel data to noaa19, using SBAFs.

    # Try I 20230404
    N19_ch1 = 0.958*M5
    N19_ch2 = 0.878*M7
    # N19_ch1 = 1.022*M5 (corr-coeff 0.9995)
    # N19_ch2 = 0.729*M6 (corr-coeff 0.9987)
    N19_ch3b = -1.1603*10^2 + 1.8362 * M12 – 1.501*10^3*M12^2 (corr-coeff 0.9942)
    N19_ch4 = 1.0003*M15 (corr-coeff 1.0!)
    N19_ch5 = 0.8 + 0.9956*M16 (corr-coeff 1.0!)

    # Try II 20231120
    r06(AVHRR,%) = 0.9393*M5 + 0.57
    r09(AVHRR,%) = 0.001275*M7*M7 + 0.7781*M7 +1.77
    tb37(AVHRR,K) = 0.9563*M12 + 9.86
    tb11(AVHRR,K) = 1.003*M15 - 0.77
    tb12(AVHRR,K) = 1.002*M16 - 0.69
    """
    
    scene["M05"].values = 0.958 * scene["M05"]
    scene["M07"].values = 0.7781 * scene["M07"] + 1.77 + 0.001275 * scene["M07"] * scene["M07"]
    scene["M15"].values = 1.003 * scene["M15"] - 0.77
    scene["M16"].values = 1.002 * scene["M16"] - 0.69
    scene["M12"].values = -1.1603 * 100 + 1.8362 * scene["M12"] - 1.501 * 0.001 * scene["M12"] * scene["M12"]

    if "npp" in scene.attrs["platform"].lower():
        scene.attrs["platform"] = "vgacsnpp"
    scene.attrs["platform"] = scene.attrs["platform"].replace("noaa", "vgac")


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
    scene.attrs['source'] = "vgac2pps.py"
    for band in REFL_BANDS:
        if band not in scene:
            continue
        # For VIIRS data sun_zenith_angle_correction_applied is applied always!
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'True'
    return nimg

def process_one_scene(scene_files, out_path, engine='h5netcdf',
                      all_channels=False, pps_channels=False, orbit_n=0, as_noaa19=False):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(
        reader='viirs_vgac_l1c_nc',
        filenames=scene_files)

    MY_MBAND = MBAND_DEFAULT

    if all_channels:
        MY_MBAND = MBANDS
    if pps_channels:
        MY_MBAND = MBAND_PPS
    if as_noaa19:
        MY_MBAND = MBAND_AVHRR

    scn_.load(MY_MBAND
              + ANGLE_NAMES
              #+ ['M12_LUT', 'M13_LUT', 'M15_LUT', 'M16_LUT']
              + ['latitude', 'longitude', 'scanline_timestamps'])

    # one ir channel
    irch = scn_['M15']

    # Set header and band attributes
    set_header_and_band_attrs(scn_, orbit_n=orbit_n)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, delete_azimuth=False)
    update_angle_attributes(scn_, irch)

    # Adjust to noaa19 with sbafs from KG
    sensor = "viirs"
    if as_noaa19:
        sensor = "avhrr"
        convert_to_noaa19(scn_)

    filename = compose_filename(scn_, out_path, instrument=sensor, band=irch)
    encoding = get_encoding_viirs(scn_)

    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor=sensor),
                       engine=engine,
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=encoding)
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
    return filename
