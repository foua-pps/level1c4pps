#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 level1c4pps developers
#
# This file is part of level1c4pps.
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>
#   Nina Hakansson <nina.hakansson@smhi.se>

"""Package Initializer for level1c4pps."""

import numpy as np
import xarray as xr
import datetime

from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

def make_azidiff_angle(sata, suna):
    """Calculate azimuth difference angle."""
    daz = abs(sata - suna)
    daz = daz % 360
    if isinstance(daz, np.ndarray):
        daz[daz > 180] = 360 - daz[daz > 180]
        return daz
    elif isinstance(daz, xr.DataArray):
        return daz.where(daz < 180, 360 - daz)
    else:
        raise ValueError("Azimuth difference is neither a Numpy nor an Xarray object! Type = %s", type(daz))

def dt64_to_datetime(dt64):
    """Conversion of numpy.datetime64 to datetime objects."""
    # https://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64/46921593#46921593
    if type(dt64) == np.datetime64:
        unix_epoch = np.datetime64(0, 's')
        one_second = np.timedelta64(1, 's')
        seconds_since_epoch = (dt64 - unix_epoch) / one_second
        dt = datetime.datetime.utcfromtimestamp(seconds_since_epoch)
        return dt
    return dt64
