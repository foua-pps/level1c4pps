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

"""Utilities to convert AVHRR GAC formattet data to PPS level-1c format."""


import os
import time
import xarray as xr
import dask.array as da
import numpy as np
from datetime import datetime
from satpy.scene import Scene
import pygac  # testing that pygac is available # noqa: F401
from level1c4pps import (get_encoding, compose_filename,
                         set_header_and_band_attrs_defaults,
                         rename_latitude_longitude, update_angle_attributes,
                         get_header_attrs, convert_angles)
import logging

logger = logging.getLogger('gac2pps')

BANDNAMES = ['1', '2', '3', '3a', '3b', '4', '5']

REFL_BANDS = ['1', '2', '3a']

PPS_TAGNAMES = {"1": "ch_r06",
                "2": "ch_r09",
                "3a": "ch_r16",
                "3": "ch_tb37",
                "3b": "ch_tb37",
                "4": "ch_tb11",
                "5": "ch_tb12"}
INSTRUMENTS = {'tirosn': 'avhrr',
               'noaa6': 'avhrr',
               'noaa7': 'avhrr/2',
               'noaa8': 'avhrr',
               'noaa9': 'avhrr/2',
               'noaa10': 'avhrr',
               'noaa11': 'avhrr/2',
               'noaa12': 'avhrr/2',
               'noaa14': 'avhrr/2',
               'noaa15': 'avhrr/3',
               'noaa16': 'avhrr/3',
               'noaa17': 'avhrr/3',
               'noaa18': 'avhrr/3',
               'noaa19': 'avhrr/3'}


def get_encoding_gac(scene):
    """Get netcdf encoding for all datasets."""
    return get_encoding(scene,
                        BANDNAMES,
                        PPS_TAGNAMES,
                        chunks=None)


def update_ancilliary_datasets(scene):
    """Rename, delete and add some datasets and attributes."""
    irch = scene['4']

    # Create new data set scanline timestamps
    first_jan_1970 = np.array([datetime(1970, 1, 1, 0, 0, 0)]).astype('datetime64[ns]')
    scanline_timestamps = np.array(scene['qual_flags'].coords['acq_time'] -
                                   first_jan_1970).astype(dtype='timedelta64[ms]').astype(np.float64)
    scene['scanline_timestamps'] = xr.DataArray(da.from_array(scanline_timestamps, chunks=1024),
                                                dims=['y'], coords={'y': scene['qual_flags']['y']})
    scene['scanline_timestamps'].attrs['units'] = 'Milliseconds since 1970-01-01 00:00:00 UTC'

    # Update qual_flags attrs
    scene['qual_flags'].attrs['id_tag'] = 'qual_flags'
    scene['qual_flags'].attrs['long_name'] = 'pygac quality flags'
    scene['qual_flags'].coords['time'] = irch.attrs['start_time']
    del scene['qual_flags'].coords['acq_time']


def set_header_and_band_attrs(scene, orbit_n=99999):
    """Set and delete some attributes."""
    irch = scene['4']
    nimg = set_header_and_band_attrs_defaults(scene, BANDNAMES, PPS_TAGNAMES, REFL_BANDS, irch, orbit_n=orbit_n)
    scene.attrs['source'] = "gac2pps.py"
    scene.attrs['is_gac'] = 'True'
    for band in BANDNAMES:
        if band not in scene:
            continue
        if band in REFL_BANDS:
            # For GAC data sun_earth_distance_correction is applied always!
            # The sun_earth_distance_correction_factor is not provided by pygac <= 1.2.1 / satpy <= 0.18.1
            scene[band].attrs['sun_earth_distance_correction_applied'] = 'True'
    return nimg


def process_one_file(gac_file, out_path='.', reader_kwargs=None, engine='h5netcdf', orbit_n=99999):
    """Make level 1c files in PPS-format."""
    tic = time.time()
    if reader_kwargs is None:
        reader_kwargs = {}
    if 'tle_dir' not in reader_kwargs:
        from pygac.configuration import get_config
        conf = get_config()
        tle_dir = conf.get('tle', 'tledir', raw=True)
        tle_name = conf.get('tle', 'tlename', raw=True)
        reader_kwargs['tle_dir'] = tle_dir
        reader_kwargs['tle_name'] = tle_name

    scn_ = Scene(reader='avhrr_l1b_gaclac',
                 filenames=[gac_file], reader_kwargs=reader_kwargs)

    # Loading all at once sometimes fails with newer satpy, so start with BANDNAMES ...

    scn_.load(BANDNAMES)
    scn_.load(['latitude',
               'longitude',
               'qual_flags',
               'sensor_zenith_angle', 'solar_zenith_angle',
               'solar_azimuth_angle', 'sensor_azimuth_angle',
               'sun_sensor_azimuth_difference_angle'])

    # one ir channel
    irch = scn_['4']

    # Set header and band attributes
    set_header_and_band_attrs(scn_, orbit_n=orbit_n)

    # Rename longitude, latitude to lon, lat.
    rename_latitude_longitude(scn_)

    # Convert angles to PPS
    convert_angles(scn_)
    update_angle_attributes(scn_, irch)

    # Handle gac specific datasets qual_flags and scanline_timestamps
    update_ancilliary_datasets(scn_)

    filename = compose_filename(scn_, out_path, instrument='avhrr', band=irch)

    encoding = get_encoding_gac(scn_)
    encoding['scanline_timestamps'].pop('units')
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_, band=irch, sensor='avhrr'),
                       engine=engine,
                       flatten_attrs=True,
                       include_lonlats=False,  # Included anyway as they are datasets in scn_
                       pretty=True,
                       encoding=encoding)

    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))
    return filename
