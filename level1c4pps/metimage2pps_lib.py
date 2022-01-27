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

"""Functions to convert EPS-SG MetImage level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import os
import time
from satpy.scene import Scene
from level1c4pps import (get_encoding, compose_filename,
                         rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         set_header_and_band_attrs_defaults,
                         convert_angles,
                         adjust_lons_to_valid_range)

import logging

# from satpy.utils import debug_on
# debug_on()

# Example:
# '/home/a001865/DATA_MISC/EPSSG_TEST/Testdata/W_xx-eumetsat-darmstadt,SAT,SGA1-VII-1B-RAD_C_EUMT_20191001043852_G_D_20070912095840_20070912100343_T_N____.nc'


logger = logging.getLogger('metimage2pps')

BANDNAMES_DEFAULT = ["vii_668",
                     "vii_865",
                     "vii_1375",
                     "vii_1630",
                     "vii_3740",
                     "vii_7325",
                     "vii_8540",
                     "vii_10690",
                     "vii_12020",
                     # "vii_443",
                     # "vii_555",
                     # "vii_752",
                     # "vii_763",
                     # "vii_914",
                     # "vii_1240",
                     "vii_2250",
                     # "vii_3959",
                     # "vii_4050",
                     "vii_6725",
                     "vii_13345"]

BANDNAMES_PPS = ["vii_668",
                 "vii_865",
                 "vii_1375",
                 "vii_1630",
                 "vii_2250",
                 "vii_3740",
                 "vii_8540",
                 "vii_10690",
                 "vii_12020"]

REFL_BANDS = ["vii_668", "vii_865", "vii_1375", "vii_1630", "vii_443",
              "vii_555", "vii_752", "vii_763", "vii_914", "vii_1240",
              "vii_2250"]

ANGLE_NAMES = ['observation_zenith', 'solar_zenith',
               'observation_azimuth', 'solar_azimuth']

PPS_TAGNAMES = {"vii_668": "ch_r06",
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

BANDNAMES = list(PPS_TAGNAMES.keys())


def get_encoding_metimage(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene, orbit_n=00000):
    """Set and delete some attributes."""
    nimg = 0  # name of first dataset is image0
    # Set some header attributes:
    irch = scene['vii_10690']
    scene.attrs['source'] = "metimage2pps.py"
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    for band in REFL_BANDS:
        if band not in scene:
            continue
        print("Is this correct, it was in testdata3.")
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'True'
    return nimg


def process_one_scene(scene_files, out_path,
                      engine='h5netcdf',
                      all_channels=False, pps_channels=False,
                      orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(reader='vii_l1b_nc', filenames=scene_files)

    MY_BANDNAMES = BANDNAMES_DEFAULT
    if all_channels:
        MY_BANDNAMES = BANDNAMES
    if pps_channels:
        MY_BANDNAMES = BANDNAMES_PPS

    scn_.load(MY_BANDNAMES + ANGLE_NAMES + ['lat_pixels', 'lon_pixels'])

    # Transpose data to get scanlines as row dimension
    for key in MY_BANDNAMES + ANGLE_NAMES + ['lat_pixels', 'lon_pixels']:
        if 'num_pixels' in scn_[key].dims:
            # satpy <= 0 .26.0
            scn_[key] = scn_[key].transpose('num_lines', 'num_pixels')
        elif scn_[key].dims[0] == 'x':
            # first dim should be y
            scn_[key] = scn_[key].transpose('y', 'x')

    # one ir channel
    irch = scn_['vii_10690']

    # Set header and band attributes
    set_header_and_band_attrs(scn_, orbit_n=orbit_n)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Adjust lons to valid range:
    adjust_lons_to_valid_range(scn_)

    # Convert angles to PPS
    convert_angles(scn_, delete_azimuth=True)
    update_angle_attributes(scn_, irch)

    # Apply sunz correction
    # apply_sunz_correction(scn_, REFL_BANDS)

    filename = compose_filename(scn_, out_path, instrument='metimage', band=irch)
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='metimage'),
                       engine=engine,
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_metimage(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
