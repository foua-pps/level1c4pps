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
                         ANGLE_ATTRIBUTES, rename_latitude_longitude,
                         update_angle_attributes, get_header_attrs,
                         convert_angles)
import pyspectral  # testing that pyspectral is available # noqa: F401
import logging

# Example:

logger = logging.getLogger('viirs2pps')

BANDNAMES = ["M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08",
             "M09", "M10", "M11", "M12", "M13", "M14", "M15", "M16",
             "I01", "I02", "I03", "I04", "I05"]

MBANDS = ["M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08",
          "M09", "M10", "M11", "M12", "M13", "M14", "M15", "M16"]
             
IBANDS = ["I01", "I02", "I03", "I04", "I05"]

REFL_BANDS = ["M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08",
              "M09", "M10", "M11", "I01", "I02", "I03"]

MBAND_PPS = [ "M05", "M07", "M09", "M10", "M12", "M14", "M15", "M16"]

IBAND_PPS = ["I01", "I02", "I03", "I04", "M09", "M14", "M15", "M16"]

ANGLE_NAMES = ['satellite_zenith_angle', 'solar_zenith_angle',
               'satellite_azimuth_angle', 'solar_azimuth_angle']

PPS_TAGNAMES = {"M05": 'ch_r06',                
                "M07": 'ch_r09',
                "M09": 'ch_r13',
                "M20": 'ch_r16',
                "M12": 'ch_tb37',
                "M14": 'ch_tb85',
                "M15": 'ch_tb11',
                "M16": 'ch_tb12',
                "I01": 'ch_r06',                
                "I02": 'ch_r09',
                "I03": 'ch_r16',
                "I04": 'ch_tb37',
}


def get_encoding_viirs(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def set_header_and_band_attrs(scene):
    """Set and delete some attributes."""
    irch = scene['M15']
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch)
    scene.attrs['source'] = "viirs2pps.py"
    # Perhaps one can get the orbit number from the h5 file?
    scene.attrs['orbit_number'] = 99999
    for band in REFL_BANDS:
        if band not in scene:
            continue
        # For VIIRS data sun_zenith_angle_correction_applied is applied always!
        scene[band].attrs['sun_zenith_angle_correction_applied'] = 'True' 
    return nimg


def process_one_scene(scene_files, out_path, use_iband_res=False):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    scn_ = Scene(
        reader='viirs_sdr',
        filenames=scene_files)
    BANDNAMES = MBAND_PPS
    resolution = 742
    if use_iband_res:
        BANDNAMES = IBAND_PPS
        resolution = 371

    scn_.load(BANDNAMES + ANGLE_NAMES, resolution=resolution)
    #import pdb;pdb.set_trace()    
    if use_iband_res:
        scn_.load(['i_latitude', 'i_longitude'])
    else:
        scn_.load(['m_latitude', 'm_longitude'])
  
        
    # one ir channel
    irch = scn_['M15']

    # Set header and band attributes
    set_header_and_band_attrs(scn_)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_, delete_azimuth=True)
    update_angle_attributes(scn_, irch)

    filename = compose_filename(scn_, out_path, instrument='viirs', band=irch)
    #import pdb;pdb.set_trace()    
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='viirs'),
                       engine='netcdf4',
                       include_lonlats=False,
                       flatten_attrs=True,
                       encoding=get_encoding_viirs(scn_))
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
