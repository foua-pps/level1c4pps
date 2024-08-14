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
#   Salomom Eliasson <salomon.eliasson@smhi.se>

"""Functions to convert VGAC level-1c data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import os
import time
from satpy.scene import Scene
from level1c4pps import (get_encoding, compose_filename,
                         set_header_and_band_attrs_defaults,
                         rename_latitude_longitude,
                         dt64_to_datetime,
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
              "M09", "M10", "M11", "I01", "I02", "I03"]

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

    # Try II Nina
    r06(AVHRR,%) = 0.9034*M5 + 1.04
    r09(AVHRR,%) = 0.8821*M7 + 0.17
    tb37(AVHRR,K) = 0.997*M12 + 0.89
    tb11(AVHRR,K) = 1.001*M15 - 0.05
    tb12(AVHRR,K) = 0.9899*M16 + 2.5

    r06(AVHRR,%) = 0.914 * M5 + 0.83
    r09(AVHRR,%) = 0.9019 * M7 - 0.18
    tb37(AVHRR,K) = 0.9962 * M12 + 0.99
    tb11(AVHRR,K) = 1.002 * M15 - 0.31
    tb12(AVHRR,K) = 0.993 * M16 + 1.62

    r06(AVHRR,%) = 0.9211*M5 + 0.85
    r09(AVHRR,%) = 0.9242*M7 -0.78
    tb37(AVHRR,K) = 0.9951*M12 + 1.18
    tb11(AVHRR,K) = 1.003*M15 - 0.55
    tb12(AVHRR,K) = 0.9977*M16 + 0.39

    KG 20230111:
    r06(AVHRR, %) = 0.9129 * M5 + 0.99
    r09(AVHRR, %) = 0.9025 * M7
    tb37(AVHRR, K) = 0.9967 * M12 + 0.85
    tb11(AVHRR, K) = 1.002 * M15 - 0.4
    tb12(AVHRR, K) = 0.9934 * M16 +1.52
    """

    scene["M05"].values = 0.9129 * scene["M05"] + 0.99
    scene["M07"].values = 0.9025 * scene["M07"]
    scene["M15"].values = 1.002 * scene["M15"] - 0.4
    scene["M16"].values = 0.9934 * scene["M16"] + 1.52
    scene["M12"].values = 0.9967 * scene["M12"] + 0.85

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

def midnight_scene(scene):
    """Check if scene passes midnight."""
    start_date = scene["M05"].attrs["start_time"].strftime("%Y%m%d")
    end_date = scene["M05"].attrs["end_time"].strftime("%Y%m%d")
    if start_date == end_date:
        return False
    return True


def get_midnight_line_nr(scene):
    """Find midnight_line, start_time and new end_time."""
    start_date = scene["M05"].attrs["start_time"].strftime("%Y-%m-%d")
    end_date = scene["M05"].attrs["end_time"].strftime("%Y-%m-%d")
    start_fine_search = len(scene['scanline_timestamps']) - 1  # As default start the fine search from end of time array
    for ind in range(0, len(scene['scanline_timestamps']), 100):
        # Search from the beginning in large chunks (100) and break when we
        # pass midnight.
        dt_obj = dt64_to_datetime(scene['scanline_timestamps'].values[:][ind])
        date_linei = dt_obj.strftime("%Y-%m-%d")
        if date_linei == end_date:
            # We just passed midnight stop and search backwards for exact line.
            start_fine_search = ind
            break
    for indj in range(start_fine_search, start_fine_search - 100, -1):
        # Midnight is in one of the previous 100 lines.
        dt_obj = dt64_to_datetime(scene['scanline_timestamps'].values[:][indj])
        date_linei = dt_obj.strftime("%Y-%m-%d")
        if date_linei == start_date:
            # We just passed midnight this is the last line for previous day.
            midnight_linenr = indj
            break
    return  midnight_linenr



def set_exact_time_and_crop(scene, start_line, end_line, time_key='scanline_timestamps'):
    """Crop datasets and update start_time end_time objects."""
    if start_line is None:
        start_line = 0
    if end_line is None:
        end_line = len(scene[time_key]) - 1
    start_time_dt64 = scene[time_key].values[start_line]
    end_time_dt64 = scene[time_key].values[end_line]
    start_time = dt64_to_datetime(start_time_dt64)
    end_time = dt64_to_datetime(end_time_dt64)
    for ds in BANDNAMES + ANGLE_NAMES + ['latitude', 'longitude', 'scanline_timestamps']:
        if ds in scene and 'nscn' in scene[ds].dims:
            scene[ds] = scene[ds].isel(nscn=slice(start_line, end_line + 1))
            try:
                # Update scene attributes to get the filenames right
                scene[ds].attrs['start_time'] = start_time
                scene[ds].attrs['end_time'] = end_time
            except TypeError:
                pass
    if start_time_dt64 != scene[time_key].values[0]:
        raise ValueError
    if end_time_dt64 != scene[time_key].values[-1]:
        raise ValueError

def split_scene_at_midnight(scene):
    """Split scenes at midnight."""
    if midnight_scene(scene):
        midnight_linenr = get_midnight_line_nr(scene)
        scene1 = scene.copy()
        scene2 = scene.copy()
        set_exact_time_and_crop(scene1, None, midnight_linenr)
        set_exact_time_and_crop(scene2, midnight_linenr + 1, None)
        return [scene1, scene2]
    return [scene]


def process_one_scene(scene_files, out_path, engine='h5netcdf',
                      all_channels=False, pps_channels=False, orbit_n=0, as_noaa19=False, avhrr_channels=False,
                      split_files_at_midnight=True):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_in = Scene(
        reader='viirs_vgac_l1c_nc',
        filenames=scene_files)

    MY_MBAND = MBAND_DEFAULT

    if all_channels:
        MY_MBAND = MBANDS
    if pps_channels:
        MY_MBAND = MBAND_PPS
    if as_noaa19:
        MY_MBAND = MBAND_AVHRR
    if avhrr_channels:
        MY_MBAND = MBAND_AVHRR

    scn_in.load(MY_MBAND
                + ANGLE_NAMES
                # + ['M12_LUT', 'M13_LUT', 'M15_LUT', 'M16_LUT']
                + ['latitude', 'longitude', 'scanline_timestamps'])
    if split_files_at_midnight:
        scenes = split_scene_at_midnight(scn_in)
    else:
        scenes = [scn_in]
    filenames = []
    for scn_ in scenes:
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
        filenames.append(filename)
    return filenames
