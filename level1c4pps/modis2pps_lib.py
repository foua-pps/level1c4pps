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

"""Functions to convert MODIS level-1 data to a NWCSAF/PPS level-1c formatet netCDF/CF file."""

import os
import time
from satpy.scene import Scene
from level1c4pps import (get_encoding, compose_filename,
                         set_header_and_band_attrs_defaults,
                         rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         convert_angles,
                         apply_sunz_correction)

import logging
from satpy.utils import debug_on
debug_on()

# Example:
# MYD03.A2014272.0215.006.2014274124038.hdf
# MYD021KM.A2014272.0215.006.2014274142955.hdf


logger = logging.getLogger('modis2pps')

# Channels pps and not problematic (8.5 )
BANDNAMES_PPS = ['1', '2', '6', '26', '20', '29', '31', '32']

# Default channel selection
BANDNAMES_DEFAULT = ['1', '2', '6', '7', '20', '26', '27', '28', '29', '31', '32', '33']


REFL_BANDS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
              '12', '13hi', '13lo', '14hi', '14lo', '15', '16', '17', '18', '19', '26']

ANGLE_NAMES = ['satellite_zenith_angle', 'solar_zenith_angle',
               'satellite_azimuth_angle', 'solar_azimuth_angle']

PPS_TAGNAMES = {'1':  'ch_r06',
                '2':  'ch_r09',
                '26': 'ch_r13',
                '6':  'ch_r16',
                '20': 'ch_tb37',
                '29': 'ch_tb85',
                '31': 'ch_tb11',
                '32': 'ch_tb12',
                # Not used yet:
                '7':  'ch_r21',
                '27': 'ch_tb67',
                '28': 'ch_tb73',
                '33': 'ch_tb133',
                '3':  'ch_rxx',
                '4':  'ch_rxx',
                '5':  'ch_rxx',
                '8':  'ch_rxx',
                '9':  'ch_rxx',
                '10': 'ch_rxx',
                '11': 'ch_rxx',
                '12': 'ch_rxx',
                '13hi': 'ch_rxx',
                '13lo': 'ch_rxx',
                '14hi': 'ch_rxx',
                '14lo': 'ch_rxx',
                '15': 'ch_rxx',
                '16': 'ch_rxx',
                '17': 'ch_rxx',
                '18': 'ch_rxx',
                '19': 'ch_rxx',
                '21': 'ch_tbxx',
                '22': 'ch_tbxx',
                '23': 'ch_tbxx',
                '24': 'ch_tbxx',
                '25': 'ch_tbxx',
                '30': 'ch_tbxx',
                '34': 'ch_tbxx',
                '35': 'ch_tbxx',
                '36': 'ch_tbxx'}

BANDNAMES = list(PPS_TAGNAMES.keys())


def get_encoding_modis(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene, orbit_n=0):
    """Set and delete some attributes."""
    irch = scene['31']
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "modis2pps.py"
    return nimg


def process_one_scene(scene_files, out_path, engine='h5netcdf', all_channels=False, pps_channels=False, orbit_n=0):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(
        reader='modis_l1b',
        filenames=scene_files)

    MY_BANDNAMES = BANDNAMES_DEFAULT
    if all_channels:
        MY_BANDNAMES = BANDNAMES
    if pps_channels:
        MY_BANDNAMES = BANDNAMES_PPS

    scn_.load(MY_BANDNAMES + ['latitude', 'longitude'] + ANGLE_NAMES, resolution=1000)
    # one ir channel
    irch = scn_['31']

    # Set header and band attributes
    set_header_and_band_attrs(scn_, orbit_n=orbit_n)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, delete_azimuth=True)
    update_angle_attributes(scn_, irch)

    # Apply sunz correction
    apply_sunz_correction(scn_, REFL_BANDS)

    filename = compose_filename(scn_, out_path, instrument='modis', band=irch)
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='modis'),
                       engine=engine,
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_modis(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
    return filename
