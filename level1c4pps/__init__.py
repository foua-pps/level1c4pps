#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Package Initializer for level1c4pps."""

import numpy as np
import xarray as xr


def make_azidiff_angle(sata, suna, fill=None):
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
