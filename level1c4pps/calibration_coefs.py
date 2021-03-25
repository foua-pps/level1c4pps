#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 level1c4pps developers
#
# This file is part of level1c4pps.py
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
# -*- coding: utf-8 -*-
# Author(s):

#  Nina.Hakansson

"""Module with calibration coefficients for SEVIRI."""

import datetime

CALIB_MODE = 'Nominal'
COEFS_MEIRINK = dict(
    MSG1=dict(
        VIS006=dict(b=24.346, a=0.3739E-3),
        VIS008=dict(b=30.989, a=0.3111E-3),
        IR_016=dict(b=22.869, a=0.0065E-3)
    ),
    MSG2=dict(
        VIS006=dict(b=21.026, a=0.2556E-3),
        VIS008=dict(b=26.875, a=0.1835E-3),
        IR_016=dict(b=21.394, a=0.0498E-3)
    ),
    MSG3=dict(
        VIS006=dict(b=19.829, a=0.5856E-3),
        VIS008=dict(b=25.284, a=0.6787E-3),
        IR_016=dict(b=23.066, a=-0.0286E-3)
    ),
    MSG4=dict(
        VIS006=dict(b=21.040, a=0.2877E-3),
        VIS008=dict(b=24.966, a=0.6074E-3),
        IR_016=dict(b=21.236, a=0.1420E-3)
    )
)

REF_DATE = datetime.date(2000, 1, 1)
REF_TIME = datetime.datetime(2000, 1, 1, 0, 0)


def calib_meirink(platform, channel, time):
    """Get MODIS-intercalibrated gain and offset for SEVIRI VIS channels.

    Reference: http://msgcpp.knmi.nl/mediawiki/index.php/MSG-SEVIRI_solar_channel_calibration

    :returns: gain, offset [mW m-2 sr-1 (cm-1)-1]
    """
    if time < REF_TIME:
        raise ValueError('Given time ({0}) is < reference time ({1})'.format(
            time, REF_TIME))

    a = COEFS_MEIRINK[platform][channel]['a']
    b = COEFS_MEIRINK[platform][channel]['b']
    delta_days = (time - REF_TIME).total_seconds() / 3600.0 / 24.0
    gain = (b + a * delta_days) / 1000.0  # micro Watts -> milli Watts
    offset = -51.0 * gain  # Space count is 51

    return gain, offset


def calib_meirink_date(platform, channel, date):
    """Get MODIS-intercalibrated gain and offset for SEVIRI VIS channels.

    Reference: http://msgcpp.knmi.nl/mediawiki/index.php/MSG-SEVIRI_solar_channel_calibration

    :returns: gain, offset [mW m-2 sr-1 (cm-1)-1]
    """
    if date < REF_DATE:
        raise ValueError('Given date ({0}) is < reference date ({1})'.format(
            date, REF_DATE))

    a = COEFS_MEIRINK[platform][channel]['a']
    b = COEFS_MEIRINK[platform][channel]['b']
    gain = (b + a*(date - REF_DATE).days) / 1000.0  # micro Watts -> milli Watts
    offset = -51.0 * gain  # Space count is 51

    return gain, offset


def get_calibration_for_time(platform, time):
    """Get MODIS-intercalibrated gain and offset for specific time."""
    coefs = {}
    for channel in ('VIS006', 'VIS008', 'IR_016'):
        gain, offset = calib_meirink(platform=platform, channel=channel,
                                     time=time)
        coefs[channel] = {'gain': gain, 'offset': offset}

    return coefs


def get_calibration_for_date(platform, date):
    """Get MODIS-intercalibrated gain and offset for specific date."""
    coefs = {}
    for channel in ('VIS006', 'VIS008', 'IR_016'):
        gain, offset = calib_meirink_date(platform=platform, channel=channel,
                                          date=date)
        coefs[channel] = {'gain': gain, 'offset': offset}

    return coefs


if __name__ == '__main__':
    time = datetime.datetime(2018, 1, 18, 12, 0)
    platform = 'MSG3'

    coefs = {}
    for channel in ('VIS006', 'VIS008', 'IR_016'):
        gain, offset = calib_meirink(platform=platform, channel=channel,
                                     time=time)
        coefs[channel] = {'gain': gain, 'offset': offset}

    print(coefs)
