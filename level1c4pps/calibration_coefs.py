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

REF_TIME = datetime.datetime(2000, 1, 1, 0, 0)
SPACE_COUNT = -51.0


def calib_meirink(platform, channel, time):
    """Get MODIS-intercalibrated gain and offset for SEVIRI VIS channels.

    Reference: http://msgcpp.knmi.nl/mediawiki/index.php/MSG-SEVIRI_solar_channel_calibration

    :returns: gain, offset [mW m-2 sr-1 (cm-1)-1]
    """
    time = _convert_to_datetime(time)
    _check_time(time)
    a = COEFS_MEIRINK[platform][channel]['a']
    b = COEFS_MEIRINK[platform][channel]['b']
    days_since_ref_time = _get_days_since_ref_time(time)
    return _calc_gain_offset(a, b, days_since_ref_time)


def _check_time(time):
    if time < REF_TIME:
        raise ValueError('Given time ({0}) is < reference time ({1})'.format(
            time, REF_TIME))


def _convert_to_datetime(date_or_time):
    if isinstance(date_or_time, datetime.date):
        return datetime.datetime.combine(date_or_time, datetime.time(0))
    return date_or_time


def _get_days_since_ref_time(time):
    return (time - REF_TIME).total_seconds() / 3600.0 / 24.0


def _calc_gain_offset(a, b, days_since_ref_time):
    gain = (b + a * days_since_ref_time)
    gain = _microwatts_to_milliwatts(gain)
    offset = SPACE_COUNT * gain
    return gain, offset


def _microwatts_to_milliwatts(microwatts):
    return microwatts / 1000.0


def get_calibration(platform, time):
    """Get MODIS-intercalibrated gain and offset for specific time."""
    coefs = {}
    for channel in ('VIS006', 'VIS008', 'IR_016'):
        gain, offset = calib_meirink(platform=platform, channel=channel,
                                     time=time)
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
