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
from sbafs_ann.convert_vgac import convert_to_vgac_with_nn
from level1c4pps import (get_encoding, compose_filename,
                         set_header_and_band_attrs_defaults,
                         rename_latitude_longitude,
                         dt64_to_datetime,
                         update_angle_attributes, get_header_attrs,
                         convert_angles)
import pyspectral  # testing that pyspectral is available # noqa: F401
import logging
import numpy as np

logger = logging.getLogger("vgac2pps")

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

# "M10", "M14" are not AVHRR channels, but needed for NN SABAF
MBAND_AVHRR = ["M05", "M07", "M12", "M15", "M16", "M10", "M14"]

MBAND_DEFAULT = ["M05", "M07", "M09", "M10",  "M11", "M12", "M14", "M15", "M16"]


ANGLE_NAMES = ["vza", "sza", "azn", "azi"]

PPS_TAGNAMES = {"M05": "ch_r06",
                "M07": "ch_r09",
                "M09": "ch_r13",
                "M10": "ch_r16",
                "M12": "ch_tb37",
                "M11": "ch_r22",
                "M14": "ch_tb85",
                "M15": "ch_tb11",
                "M16": "ch_tb12",
                "I01": "ch_r06",
                "I02": "ch_r09",
                "I03": "ch_r16",
                "I04": "ch_tb37",
                # Not used by pps:
                "I05": "ch_tbxx",
                "M01": "ch_rxx",
                "M02": "ch_rxx",
                "M03": "ch_rxx",
                "M04": "ch_rxx",
                "M06": "ch_rxx",
                "M08": "ch_rxx",
                "M13": "ch_tbxx"}

# SBAF dictionary
#
# Spectral band adjustment factors for converting VGAC to AVHRR

SBAF = {
    "v2": {
        "r06": {
            "viirs_channel": "M05",
            "slope": 0.9129,
            "offset": 0.99,
            "comment": "based on nadir collocation data for SZA < 75",
        },
        "r09": {
            "viirs_channel": "M07",
            "slope": 0.9025,
            "offset": 0,
            "comment": "based on nadir collocation data for SZA < 75",
        },
        "tb37": {
            "viirs_channel": "M12",
            "slope": 0.9967,
            "offset": 0.85,
            "comment": "based on nadir collocation data for SZA 95-180",
        },
        "tb11": {
            "viirs_channel": "M15",
            "slope": 1.002,
            "offset": -0.4,
            "comment": "based on nadir collocation data for SZA 0 -180",
        },
        "tb12": {
            "viirs_channel": "M16",
            "slope": 0.9934,
            "offset":  1.52,
            "comment": "based on nadir collocation data for SZA 0-180",
        },
    },
    "v3": {
        "r06": {
            "viirs_channel": "M05",
            "slope": 0.8949,
            "offset": 2.07,
            "comment": "based on nadir collocation data for SZA < 75",
        },
        "r09": {
            "viirs_channel": "M07",
            "slope": 0.9204,
            "offset": -0.87,
            "comment": "based on nadir collocation data for SZA < 75",
        },
        "tb37_night": {
            "viirs_channel": "M12",
            "slope": 0.9923,
            "offset": 1.8,
            "min_sunzenith": 89,
            "comment": "based on nadir collocation data for SZA 100-180. SZA >= 89",
        },
        "tb37_day": {
            "viirs_channel": "M12",
            "slope": 0.9491,
            "offset": 16.0,
            "max_sunzenith": 80,
            "comment": "based on NASA SBAFS for SZA < 75 o and VZA < 5. SZA < 80",
        },
        "tb37_twilight": {
            "viirs_channel": "M12",
            "slope": 0.9707,
            "offset": 8.9,
            "min_sunzenith": 80,
            "max_sunzenith": 89,
            "comment": "The linear average of the SBAFs for t37_day and t37_night. 80<= SZA <89",
        },
        "tb11": {
            "viirs_channel": "M15",
            "slope": 1.003,
            "offset": -0.78,
            "comment": "based on nadir collocation data for SZA 0 -180",
        },
        "tb12": {
            "viirs_channel": "M16",
            "slope": 1.003,
            "offset": -0.84,
            "comment": "based on nadir collocation data for SZA 0-180",
            },
    },
    "v4": {
        "r06": {
            "viirs_channel": "M05",
            "slope": 0.8949,
            "offset": 2.07,
            "comment": "based on nadir collocation data for SZA < 75",
        },
        "r09": {
            "viirs_channel": "M07",
            "slope": 0.9204,
            "offset": -0.87,
            "comment": "based on nadir collocation data for SZA < 75",
        },
        "tb37_night": {
            "viirs_channel": "M12",
            "slope": 0.9923,
            "offset": 1.8,
            "min_sunzenith": 89,
            "comment": "based on nadir collocation data for SZA 100-180. SZA >= 89",
        },
        "tb37_day": {
            "viirs_channel": "M12",
            "slope": 0.9491,
            "offset": 16.0,
            "max_sunzenith": 80,
            "comment": "based on NASA SBAFS for SZA < 75 o and VZA < 5. SZA < 80",
        },
        "tb37_twilight": {
            "viirs_channel": "M12",
            "slope": 0.9707,
            "offset": 8.9,
            "min_sunzenith": 80,
            "max_sunzenith": 89,
            "comment": "The linear average of the SBAFs for t37_day and t37_night. 80<= SZA <89",
        },
        "tb11": {
            "viirs_channel": "M15",
            "slope": 1.002,
            "offset": -0.4,
            "comment": "based on nadir collocation data for SZA 0 -180",
        },
        "tb12": {
            "viirs_channel": "M16",
            "slope": 0.9934,
            "offset": 1.52,
            "comment": "based on nadir collocation data for SZA 0-180",
        },
    },
    "v5": {
        "r06": {
            "viirs_channel": "M05",
            "slope": 0.8949,
            "offset": 2.07,
            "comment": "based on nadir collocation data for SZA < 75",
        },
        "r09": {
            "viirs_channel": "M07",
            "slope": 0.9204,
            "offset": -0.87,
            "comment": "based on nadir collocation data for SZA < 75",
        },
        "tb37_night": {
            "viirs_channel": "M12",
            "slope": 0.9967,
            "offset": 0.85,
            "min_sunzenith": 89,
            "comment": "based on nadir collocation data for SZA >= 89",
        },
        "tb37_day": {
            "viirs_channel": "M12",
            "slope": 0.9455,
            "offset": 12.69,
            "max_sunzenith": 80,
            "comment": "based on NASA SBAFS for VZA < 5 and SZA < 80",
        },
        "tb37_twilight": {
            "viirs_channel": "M12",
            "slope": 0.9711,
            "offset": 6.77,
            "min_sunzenith": 80,
            "max_sunzenith": 89,
            "comment": "The linear average of the SBAFs for t37_day and t37_night. 80<= SZA <89",
        },
        "tb11": {
            "viirs_channel": "M15",
            "slope": 1.002,
            "offset": -0.4,
            "comment": "based on nadir collocation data for SZA 0 -180",
        },
        "tb12": {
            "viirs_channel": "M16",
            "slope": 0.9934,
            "offset": 1.52,
            "comment": "based on nadir collocation data for SZA 0-180",
        },
    },
    "v6": {
        "r06": {
            "viirs_channel": "M05",
            "slope": 0.9393,
            "offset": 0.57,
            "comment": "Based on collocation data, VZA < 3 for SZA < 60",
        },
        "r09": {
            "viirs_channel": "M07",
            "slope": 0.9104,
            "offset": -0.90,
            "comment": "Based on collocation data, VZA < 3 for SZA < 60",
        },
        "tb37_night": {
            "viirs_channel": "M12",
            "slope": 0.9967,
            "offset": 0.85,
            "min_sunzenith": 89,
            "comment": "Based on collocations data for SZA >= 89",
        },
        "tb37_day": {
            "viirs_channel": "M12",
            "slope": 0.9455,
            "offset": 12.69,
            "max_sunzenith": 80,
            "comment": "Based on collocations for SZA < 80 o and VZA < 5",
        },
        "tb37_twilight": {
            "viirs_channel": "M12",
            "slope": 0.9711,
            "offset": 6.77,
            "min_sunzenith": 80,
            "max_sunzenith": 89,
            "comment": "The linear average of the SBAFs for t37_day and t37_night. 80<= SZA <89",
        },
        "tb11": {
            "viirs_channel": "M15",
            "slope": 1.003,
            "offset": -0.77,
            "comment": "Based on collocation data for VZA < 3 and SZA 0 -180",
        },
        "tb12": {
            "viirs_channel": "M16",
            "slope": 1.002,
            "offset": -0.69,
            "comment": "Based on collocation data for VZA < 3 and SZA 0-180",
        }
    },
    "v7": {
        "r06": {
            "viirs_channel": "M05",
            "slope": 0.8534,
            "offset":1.8517,
            "comment": "Based on collocation data for VZA < 15",
        },
        "r09": {
            "viirs_channel": "M07",
            "slope":0.8507,
            "offset":1.1157,
            "comment": "Based on collocation data for VZA < 15",
        },
        "tb37": {
            "viirs_channel": "M12",
            "slope":0.9734,
            "offset":6.1707,
            "comment": "Based on collocation data for VZA < 15",
        },
        "tb11": {
            "viirs_channel": "M15",
            "slope": 1.0006,
            "offset": -0.0378,
            "comment": "Based on collocation data for VZA < 15",
        },
        "tb12": {
            "viirs_channel": "M16",
            "slope": 0.9906,
            "offset": 2.1505,
            "comment": "Based on collocation data for VZA < 15",
        }
    },
    "v8": {
        "r06_day": {
            "viirs_channel": "M05",
            "slope": 0.8960,
            "offset": 1.7118,
            "max_sunzenith": 80,
            "comment": "Based on collocation data for VZA < 15",
        },
        "r09_day": {
            "viirs_channel": "M07",
            "slope": 0.8907,
            "offset": 0.5547,
            "max_sunzenith": 80,
            "comment": "Based on collocation data for VZA < 15",
        },
        "tb37_day": {
            "viirs_channel": "M12",
            "slope": 0.9817,
            "offset": 3.0744,
            "max_sunzenith": 80,
            "comment": "Based on collocations for VZA < 15 and SZA < 80",
        },
        "tb11_day": {
            "viirs_channel": "M15",
            "slope": 0.9926,
            "offset": 2.0934,
            "max_sunzenith": 80,
            "comment": "Based on collocation data for VZA < 15 and SZA < 80",
        },
        "tb12_day": {
            "viirs_channel": "M16",
            "slope": 0.9774,
            "offset": 5.6555,
            "max_sunzenith": 80,
            "comment": "Based on collocation data for VZA < 15 and SZA < 80",
        },
        "r06_twilight": {
            "viirs_channel": "M05",
            "slope": 0.6710,
            "offset": 4.5075,
            "min_sunzenith": 80,
            "max_sunzenith": 90,
            "comment": "Based on collocation data for VZA < 15 and 80<= SZA <=90",
        },
        "r09_twilight": {
            "viirs_channel": "M07",
            "slope": 0.7120,
            "offset": 3.7473,
            "min_sunzenith": 80,
            "max_sunzenith": 90,
            "comment": "Based on collocation data for VZA < 15 and 80<= SZA <=90",
        },
        "tb37_twilight": {
            "viirs_channel": "M12",
            "slope": 0.9711,
            "offset": 6.77,
            "min_sunzenith": 80,
            "max_sunzenith": 90,
            "comment": "Based on collocation data for VZA < 15 and 80<= SZA <=90",
        },
        "tb11_twilight": {
            "viirs_channel": "M15",
            "slope": 1.0024,
            "offset": -0.5164,
            "min_sunzenith": 80,
            "max_sunzenith": 90,
            "comment": "Based on collocation data for VZA < 15 and 80<= SZA <=90",
        },
        "tb12_twilight": {
            "viirs_channel": "M16",
            "slope": 0.9973,
            "offset": -0.7479,
            "min_sunzenith": 80,
            "max_sunzenith": 90,
            "comment": "Based on collocation data for VZA < 15 and 80<= SZA <=90",
        },
        "tb37_night": {
            "viirs_channel": "M12",
            "slope": 0.9948,
            "offset": 1.3254,
            "min_sunzenith": 90,
            "comment": "Based on collocations for VZA < 15  and SZA >90",
        },
        "tb11_night": {
            "viirs_channel": "M15",
            "slope": 1.0042,
            "offset": -0.9186,
            "min_sunzenith": 90,
            "comment": "Based on collocations for VZA < 15  and SZA >90",
        },
        "tb12_night": {
            "viirs_channel": "M16",
            "slope": 0.9962,
            "offset": 0.7373,
            "min_sunzenith": 90,
            "comment": "Based on collocation data for VZA < 15 and SZA > 90",
        },
    },
    "KNMI": {
        "r06": {
            "viirs_channel": "M05",
            "slope": 0.9401,
            "offset": 0.629,
            "comment": "VZA<60, SZA<60, Delta(VZA) < 5, Delta(Scat-angle) < 5",
        },
        "r09": {
            "viirs_channel": "M07",
            "slope": 0.9345,
            "offset": -1.304,
            "comment": "VZA<60, SZA<60, Delta(VZA) < 5, Delta(Scat-angle) < 5",
        },
        "tb37_night": {
            "viirs_channel": "M12",
            "slope": 0.9934,
            "offset": 1.659,
            "min_sunzenith": 89,
            "comment": " VZA<60, SZA>95, Delta(VZA) < 5",
        },
        "tb37_day": {
            "viirs_channel": "M12",
            "slope": 0.9572,
            "offset": 9.572,
            "max_sunzenith": 80,
            "comment": "VZA<60, SZA<60, Delta(VZA) < 5, Delta(Scat-angle) < 5",
        },
        "tb37_twilight": {
            "viirs_channel": "M12",
            "slope": 0.9753,
            "offset": 5.6155,
            "min_sunzenith": 80,
            "max_sunzenith": 89,
            "comment": "The linear average of the SBAFs for t37_day and t37_night. 80<= SZA <89",
        },
        "tb11": {
            "viirs_channel": "M15",
            "slope": 0.9986,
            "offset": 0.600,
            "comment": "",
        },
        "tb12": {
            "viirs_channel": "M16",
            "slope": 0.9943,
            "offset": 1.328,
            "comment": "",
        },
    },
    "KNMI_v2": {
        "r06": {
            "viirs_channel": "M05",
            "slope": 0.9401,
            "offset": 0.629,
            "comment": "VZA<60, SZA<60, Delta(VZA) < 5, Delta(Scat-angle) < 5",
        },
        "r09": {
            "viirs_channel": "M07",
            "slope": 0.9345,
            "offset": -1.304,
            "comment": "VZA<60, SZA<60, Delta(VZA) < 5, Delta(Scat-angle) < 5",
        },
        "tb37_day": {
            "viirs_channel": "M12",
            "slope": 0.9572,
            "offset": 9.572,
            "max_sunzenith": 70,
            "comment": "VZA<60, SZA<60, Delta(VZA) < 5, Delta(Scat-angle) < 5",
        },
        "tb37_night": {
            "viirs_channel": "M12",
            "slope": 0.9934,
            "offset": 1.659,
            "min_sunzenith": 85,
            "comment": "VZA<60, SZA>95, Delta(VZA) < 5",
        },
        "tb37_twilight": {
            "viirs_channel": "M12",
            "slope": None,
            "offset": None,
            "min_sunzenith": 70,
            "max_sunzenith": 85,
            "comment": "70 < SZA < 85: BT = (1-f)*BT(day) + f*BT(night), f=(SZA-70)/15",
        },
        "tb11": {
            "viirs_channel": "M15",
            "slope": 0.9986,
            "offset": 0.600,
            "comment": "VZA<60, SZA>95, Delta(VZA) < 5",
        },
        "tb12": {
            "viirs_channel": "M16",
            "slope": None,
            "offset": None,
            "comment": "VZA<60, SZA>95, Delta(VZA) < 5. Applies SBAF(tb11)-1.1646*(tb11-tb12)-0.235",
        },
    },
    "NN_v1": {
        "cfg_file_day": "ch7_SATZ_less_15_SUNZ_0_89_TD_1_min.yaml",
        "cfg_file_night": "ch4_SATZ_less_25_SUNZ_90_180_TD_5_min.yaml",
        "cfg_file_twilight": None,
        "comment": "NN based on AVHRR and VGAC matchups using all AVHRR heritage channels"
        },
    "NN_v2": {
        "cfg_file_day": "ch7_satz_max_25_SUNZ_0_80_tdiff_300_sec_20241031.yaml",
        "cfg_file_night": "ch4_satz_max_25_SUNZ_90_180_tdiff_300_sec_20241031.yaml",
        "cfg_file_twilight": "ch7_satz_max_25_SUNZ_80_89_tdiff_300_sec_20241031.yaml",
        "comment": "NN based on AVHRR and VGAC matchups using all AVHRR heritage channels"
        },
    "NN_v3": {
        "cfg_file_day": "ch7_satz_max_15_SUNZ_0_80_tdiff_120_sec_20241120.yaml",
        "cfg_file_night": "ch4_satz_max_15_SUNZ_90_180_tdiff_120_sec_20241120.yaml",
        "cfg_file_twilight": "ch7_satz_max_15_SUNZ_80_89_tdiff_120_sec_20241120.yaml",
        "comment": "NN based on AVHRR and VGAC matchups using all AVHRR heritage channels"
        },
    "NN_v4": {
        "cfg_file_day": "ch7_satz_max_15_SUNZ_0_80_tdiff_120_sec_20241204.yaml",
        "cfg_file_night": "ch4_satz_max_15_SUNZ_90_180_tdiff_120_sec_20241204.yaml",
        "cfg_file_twilight": "ch7_satz_max_15_SUNZ_80_89_tdiff_120_sec_20241204.yaml",
        "comment": "NN based on AVHRR and VGAC matchups using all AVHRR heritage channels"
    }
}


def convert_to_noaa19_neural_network(scene, sbaf_version):
    """Applies AVHRR SBAF to VGAC channels using NN approach"""

    if sbaf_version in ["NN_v1", "NN_v2", "NN_v3"]:
        day_cfg_file = SBAF[sbaf_version]['cfg_file_day']
        night_cfg_file = SBAF[sbaf_version]['cfg_file_night']
        twilight_cfg_file = SBAF[sbaf_version]['cfg_file_twilight']
    else:
        logger.exception(f"Unrecognized NN version, {sbaf_version}")
    scene = convert_to_vgac_with_nn(scene, day_cfg_file, night_cfg_file, twilight_cfg_file)

    logger.info(f'Created NN version {sbaf_version}')


def convert_to_noaa19_linear(scene, SBAF):
    """ Apply linear regression"""

    for avhhr_chan, scaling in SBAF.items():
        viirs_channel = scaling["viirs_channel"]
        offset = scaling["offset"]
        comment = scaling["comment"]
        slope = scaling["slope"]
        filt = np.ones_like(scene[viirs_channel].values, dtype=bool)
        if "min_sunzenith" in scaling:
            filt = filt & (scene["sunzenith"].values >= scaling["min_sunzenith"])
        if "max_sunzenith" in scaling:
            filt = filt & (scene["sunzenith"].values < scaling["max_sunzenith"])
        scene[viirs_channel].values = np.where(filt,
                                               slope * scene[viirs_channel].values + offset,
                                               scene[viirs_channel].values
                                               )
        logger.info(f"{avhhr_chan:<13} = {slope:<6}*{viirs_channel:<3}+{offset:<5} ({comment})")


def convert_to_noaa19_KNMI_v2(scene, sbaf_version):
    """ Apply 1 channel linear regression SBAF for KNMI version 2"""

    # I need to save the t11 values before the SBAF adjustment as they are needed for the tb12 SBAF
    tb11_original = scene["M15"].values.copy()
    for avhrr_chan, scaling in SBAF[sbaf_version].items():
        viirs_channel = scaling["viirs_channel"]
        offset = scaling["offset"]
        comment = scaling["comment"]
        slope = scaling["slope"]
        filt = np.ones_like(scene[viirs_channel].values, dtype=bool)
        if "min_sunzenith" in scaling:
            filt = filt & (scene["sunzenith"].values >= scaling["min_sunzenith"])
        if "max_sunzenith" in scaling:
            filt = filt & (scene["sunzenith"].values < scaling["max_sunzenith"])

        if slope is not None:
            # Then simple linear regression
            scene[viirs_channel].values = np.where(
                filt,
                slope * scene[viirs_channel].values + offset,
                scene[viirs_channel].values
                )
            logger.info(f"{avhrr_chan:<13} = {slope:<6}*{viirs_channel:<3}+{offset:<5} ({comment})")
        else:
            if avhrr_chan == "tb37_twilight":
                # 70 < SZA < 85: BT = (1-f)*BT(day) + f*BT(night), f=(SZA-70)/15

                f = (scene["sunzenith"].values - scaling["min_sunzenith"])/15
                tb37_day_slope = scaling["tb37_day"]["slope"]
                tb37_day_offset = scaling["tb37_day"]["offset"]
                tb37_night_slope = scaling["tb37_night"]["slope"]
                tb37_night_offset = scaling["tb37_night"]["offset"]
                tb37_day = tb37_day_slope * scene[viirs_channel].values + tb37_day_offset
                tb37_night = tb37_night_slope * scene[viirs_channel].values + tb37_night_offset

                scene[viirs_channel].values = np.where(
                    filt,
                    (1-f)*tb37_day + f*tb37_night,
                    scene[viirs_channel].values
                )

            elif avhrr_chan == "tb12":
                # BT(5) = BT(ch4)-1.1646*(BT(M15)-BT(M16))-0.235
                scene[viirs_channel].values = scene["M15"].values-1.1646*(tb11_original-scene["M16"].values)-0.235
            else:
                logger.exception(f'Unknown channel, {avhrr_chan}, or missing slope parameter')
            logger.info(f"{avhrr_chan:<13}: ({comment})")


def convert_to_noaa19(scene, sbaf_version):
    """ Applies AVHRR SBAF to VGAC channels"""

    logger.info(f"Using SBAF_{sbaf_version}")

    if "NN" in sbaf_version:
        convert_to_noaa19_neural_network(scene, sbaf_version)
    elif sbaf_version == "KNMI_v2":
        convert_to_noaa19_KNMI_v2(scene, sbaf_version)
    else:
        convert_to_noaa19_linear(scene, SBAF[sbaf_version])

    if "npp" in scene.attrs["platform"].lower():
        scene.attrs["platform"] = "vgacsnpp"
    scene.attrs["platform"] = scene.attrs["platform"].replace("noaa", "vgac")

    # These are only needed for the NN, but they must be removed before saving
    del scene["M10"]
    del scene["M14"]


def get_encoding_viirs(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene, orbit_n=0):
    """Set and delete some attributes."""
    irch = scene["M15"]
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    scene.attrs["source"] = "vgac2pps.py"
    for band in REFL_BANDS:
        if band not in scene:
            continue
        # For VIIRS data sun_zenith_angle_correction_applied is applied always!
        scene[band].attrs["sun_zenith_angle_correction_applied"] = "True"
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
    start_fine_search = len(scene["scanline_timestamps"]) - 1  # As default start the fine search from end of time array
    file_contain_bad_time_info = False
    for ind in range(0, len(scene["scanline_timestamps"]), 100):
        # Search from the beginning in large chunks (100) and break when we
        # pass midnight.
        if np.isnan(scene["scanline_timestamps"].values[:][ind]):
            # Sometimes time info is wrong 10^36 hours since ...
            file_contain_bad_time_info = True
            continue
        dt_obj = dt64_to_datetime(scene["scanline_timestamps"].values[:][ind])
        date_linei = dt_obj.strftime("%Y-%m-%d")
        if date_linei == end_date:
            # We just passed midnight stop and search backwards for exact line.
            start_fine_search = ind
            break
    for indj in range(start_fine_search, start_fine_search - 100, -1):
        # Midnight is in one of the previous 100 lines.
        if np.isnan(scene["scanline_timestamps"].values[:][indj]):
            raise ValueError("Error in time information in VGAC file.")
        dt_obj = dt64_to_datetime(scene["scanline_timestamps"].values[:][indj])
        date_linei = dt_obj.strftime("%Y-%m-%d")
        if date_linei == start_date:
            # We just passed midnight this is the last line for previous day.
            midnight_linenr = indj
            break
        if file_contain_bad_time_info and indj == start_fine_search - 99:
            raise ValueError("Error in time information in VGAC file.")
    return midnight_linenr


def set_exact_time_and_crop(scene, start_line, end_line, time_key="scanline_timestamps"):
    """Crop datasets and update start_time end_time objects."""
    if start_line is None:
        start_line = 0
    if end_line is None:
        end_line = len(scene[time_key]) - 1
    start_time_dt64 = scene[time_key].values[start_line]
    end_time_dt64 = scene[time_key].values[end_line]
    start_time = dt64_to_datetime(start_time_dt64)
    end_time = dt64_to_datetime(end_time_dt64)
    for ds in BANDNAMES + ANGLE_NAMES + ["latitude", "longitude", "scanline_timestamps"]:
        if ds in scene and "nscn" in scene[ds].dims:
            scene[ds] = scene[ds].isel(nscn=slice(start_line, end_line + 1))
            try:
                # Update scene attributes to get the filenames right
                scene[ds].attrs["start_time"] = start_time
                scene[ds].attrs["end_time"] = end_time
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


def fix_timestamp_datatype(scene, encoding):
    """Fix time datatype."""
    param = "scanline_timestamps"
    if "milliseconds" in encoding[param]["units"]:
        scene[param].data = scene[param].data.astype('datetime64[ms]')


def process_one_scene(scene_files, out_path, engine="h5netcdf",
                      all_channels=False, pps_channels=False, orbit_n=0,
                      noaa19_sbaf_version=None, avhrr_channels=False,
                      split_files_at_midnight=True):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_in = Scene(
        reader="viirs_vgac_l1c_nc",
        filenames=scene_files)

    MY_MBAND = MBAND_DEFAULT

    if all_channels:
        MY_MBAND = MBANDS
    if pps_channels:
        MY_MBAND = MBAND_PPS
    if noaa19_sbaf_version is not None:
        MY_MBAND = MBAND_AVHRR
    if avhrr_channels:
        MY_MBAND = MBAND_AVHRR

    scn_in.load(MY_MBAND
                + ANGLE_NAMES
                + ["latitude", "longitude", "scanline_timestamps"])
    if split_files_at_midnight:
        scenes = split_scene_at_midnight(scn_in)
    else:
        scenes = [scn_in]
    filenames = []
    for scn_ in scenes:
        # one ir channel
        irch = scn_["M15"]

        # Set header and band attributes
        set_header_and_band_attrs(scn_, orbit_n=orbit_n)

        # Rename longitude, latitude to lon, lat.
        rename_latitude_longitude(scn_)

        # Convert angles to PPS
        convert_angles(scn_, delete_azimuth=False)
        update_angle_attributes(scn_, irch)
        # Adjust to noaa19 with sbafs from KG
        sensor = "viirs"
        if noaa19_sbaf_version is not None:
            sensor = "avhrr"
            convert_to_noaa19(scn_, noaa19_sbaf_version)

        filename = compose_filename(scn_, out_path, instrument=sensor, band=irch)
        encoding = get_encoding_viirs(scn_)
        fix_timestamp_datatype(scn_, encoding)

        scn_.save_datasets(writer="cf",
                           filename=filename,
                           header_attrs=get_header_attrs(scn_, band=irch, sensor=sensor,
                                                         sbaf_version=noaa19_sbaf_version),
                           engine=engine,
                           include_lonlats=False,
                           flatten_attrs=True,
                           encoding=encoding)
        logger.info("Saved file {:s} after {:3.1f} seconds".format(
             os.path.basename(filename),
             time.time()-tic))
        filenames.append(filename)
    return filenames
