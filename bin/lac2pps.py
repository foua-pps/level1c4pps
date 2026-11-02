#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Script to make AVHRR LAC and FRAC level1c in PPS-format with pytroll."""

import argparse

from level1c4pps.lac2pps_lib import process_one_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make a PPS level1c file from an AVHRR level 1b file.")
    parser.add_argument('file', type=str, help='AVHRR level 1b file to process')
    parser.add_argument('-o', '--out_dir', type=str, default='.', help="Output directory where to store level1c file.")
    parser.add_argument('-td', '--tle_dir', type=str, default='.', help="TLE directory.")
    parser.add_argument('-tn', '--tle_name', type=str, default='.', help="TLE file name.")
    options = parser.parse_args()
    process_one_file(options.file, options.out_dir,
                     reader_kwargs={'tle_dir': options.tle_dir, 'tle_name': options.tle_name})
