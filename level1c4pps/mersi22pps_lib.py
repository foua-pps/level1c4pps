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
from datetime import datetime
from satpy.scene import Scene
from level1c4pps import (make_azidiff_angle, get_encoding, compose_filename,
                         ANGLE_ATTRIBUTES, rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs)
import pyspectral  # testing that pyspectral is available # noqa: F401
import logging

# Example:
# tf2019234102243.FY3D-X_MERSI_GEOQK_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_GEO1K_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_1000M_L1B.HDF
# tf2019234102243.FY3D-X_MERSI_0250M_L1B.HDF
#

logger = logging.getLogger('mersi22pps')

PLATFORM_SHORTNAMES = {'FY3D': 'FY-3D'}
# BANDNAMES = ['%d' % (chn+1) for chn in range(25)]
BANDNAMES = ['3', '4', '5', '6', '20', '23', '24', '25']

PPS_TAGNAMES = {'3': 'ch_r06',
                '4': 'ch_r09',
                '5': 'ch_r13',
                '6': 'ch_r16',
                '20': 'ch_tb37',
                '22': 'ch_tb73',
                '23': 'ch_tb85',
                '24': 'ch_tb11',
                '25': 'ch_tb12'}

SATPY_ANGLE_NAMES = {
    'sunzenith': 'solar_zenith_angle',
    'satzenith': 'satellite_zenith_angle',
    'sunazimuth': 'solar_azimuth_angle',
    'satazimuth': 'satellite_azimuth_angle'}

MERSI2_LEVEL1_FILE_PATTERN = 'tf{start_time:%Y%j%H%M%S}.{platform_shortname:4s}-X_MERSI_{dataset}_L1B.HDF'


def get_encoding_mersi2(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def convert_angles(scene, satpy_angle_names):
    """Convert angles to pps format."""
    for angle in ['sunzenith', 'satzenith', 'sunazimuth', 'satazimuth']:
        scene[angle] = scene[satpy_angle_names[angle]]
        del scene[satpy_angle_names[angle]]
    scene['azimuthdiff'] = make_azidiff_angle(scene['satazimuth'], scene['sunazimuth'])
    scene['azimuthdiff'].attrs = scene['sunazimuth'].attrs
    del scene['satazimuth']
    del scene['sunazimuth']


def set_header_and_band_attrs(scene):
    """Set and delete some attributes."""
    nimg = 0  # name of first dataset is image0
    # Set some header attributes:
    irch = scene['24']
    scene.attrs['platform'] = irch.attrs['platform_name']
    sensor_name = [x for x in scene.attrs['sensor']][0]
    scene.attrs['instrument'] = sensor_name.upper()
    scene.attrs['source'] = "mersi22pps.py"
    # Perhaps one can get the orbit number from the hdf file?
    # FIXME!
    scene.attrs['orbit_number'] = "99999"
    nowutc = datetime.utcnow()
    scene.attrs['date_created'] = nowutc.strftime("%Y-%m-%dT%H:%M:%SZ")

    # bands
    for band in BANDNAMES:
        idtag = PPS_TAGNAMES.get(band, None)
        if not idtag:
            continue
        scene[band].attrs['id_tag'] = idtag
        scene[band].attrs['description'] = 'MERSI-2 band ' + str(band)
        scene[band].attrs['sun_earth_distance_correction_applied'] = 'False'
        scene[band].attrs['sun_earth_distance_correction_factor'] = 1.0
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
        scene[band].attrs['name'] = "image{:d}".format(nimg)
        scene[band].attrs['coordinates'] = 'lon lat'
        # Add time coordinate. To make cfwriter aware that we want 3D data.
        scene[band].coords['time'] = irch.attrs['start_time']
        try:
            del scene[band].attrs['area']
        except KeyError:
            pass
        nimg += 1
    return nimg


def process_one_scene(scene_files, out_path):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(
        reader='mersi2_l1b',
        filenames=scene_files)

    scn_.load(BANDNAMES + ['latitude', 'longitude',
                           'satellite_zenith_angle', 'solar_zenith_angle',
                           'satellite_azimuth_angle', 'solar_azimuth_angle'], resolution=1000)
    # one ir channel
    irch = scn_['24']

    # Set header and band attributes
    set_header_and_band_attrs(scn_)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, SATPY_ANGLE_NAMES)
    update_angle_attributes(scn_, irch)
    for angle in ['sunzenith', 'satzenith', 'azimuthdiff']:
        scn_[angle].attrs['file_key'] = ANGLE_ATTRIBUTES['mersi2_file_key'][angle]

    filename = compose_filename(scn_, out_path, instrument='mersi2', band=irch)
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='mersi-2'),
                       engine='netcdf4',
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_mersi2(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
