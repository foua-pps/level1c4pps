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
from enum import Enum


class CalibrationData(Enum):
    COEFS = dict(
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
    SPACE_COUNT = -51.0
    REF_TIME = datetime.datetime(2000, 1, 1, 0, 0)
    TIME_COVERAGE = {
        "start": datetime.datetime(2004, 1, 1, 0, 0),
        "end": datetime.datetime(2021, 1, 1, 0, 0)
    }
    SATPY_CALIB_MODE = 'Nominal'


def get_calibration(platform, time, clip=True):
    """Get MODIS-intercalibrated gain and offset for specific time.

    Args:
        platform: Platform name.
        time: Date or time of observations to be calibrated.
        clip: If True, do not extrapolate calibration coefficients beyond the
            time coverage of the calibration dataset. Instead, clip at the
            boundaries, that means return the boundary coefficients for
            timestamps outside the coverage.
    """
    coefs = {}
    for channel in ('VIS006', 'VIS008', 'IR_016'):
        coefs[channel] = _get_single_channel_calibration(
            platform=platform,
            channel=channel,
            time=time,
            clip=clip
        )
    return coefs


def _get_single_channel_calibration(platform, channel, time, clip):
    time = _prepare_time(time, clip)
    gain, offset = calib_meirink(platform, channel, time)
    return {'gain': gain, 'offset': offset}


def _prepare_time(time, clip):
    time = _convert_to_datetime(time)
    _check_time(time)
    if clip:
        time = _clip_at_coverage_bounds(time)
    return time


def _convert_to_datetime(date_or_time):
    if isinstance(date_or_time, datetime.date):
        return datetime.datetime.combine(date_or_time, datetime.time(0))
    return date_or_time


def _check_time(time):
    ref_time = CalibrationData.REF_TIME.value
    if time < ref_time:
        raise ValueError('Given time ({0}) is < reference time ({1})'.format(
            time, ref_time))


def _clip_at_coverage_bounds(time):
    time_cov = CalibrationData.TIME_COVERAGE.value
    time = max(time, time_cov["start"])
    time = min(time, time_cov["end"])
    return time


def calib_meirink(platform, channel, time):
    """Get MODIS-intercalibrated gain and offset for SEVIRI VIS channels.

    Reference: https://msgcpp.knmi.nl/solar-channel-calibration.html

    :returns: gain, offset [mW m-2 sr-1 (cm-1)-1]
    """
    coefs = CalibrationData.COEFS.value
    a = coefs[platform][channel]['a']
    b = coefs[platform][channel]['b']
    days_since_ref_time = _get_days_since_ref_time(time)
    return _calc_gain_offset(a, b, days_since_ref_time)


def _get_days_since_ref_time(time):
    ref_time = CalibrationData.REF_TIME.value
    return (time - ref_time).total_seconds() / 3600.0 / 24.0


def _calc_gain_offset(a, b, days_since_ref_time):
    gain = (b + a * days_since_ref_time)
    gain = _microwatts_to_milliwatts(gain)
    offset = CalibrationData.SPACE_COUNT.value * gain
    return gain, offset


def _microwatts_to_milliwatts(microwatts):
    return microwatts / 1000.0


if __name__ == '__main__':
    time = datetime.datetime(2018, 1, 18, 12, 0)
    platform = 'MSG3'

    coefs = {}
    for channel in ('VIS006', 'VIS008', 'IR_016'):
        gain, offset = calib_meirink(platform=platform, channel=channel,
                                     time=time)
        coefs[channel] = {'gain': gain, 'offset': offset}

    print(coefs)
