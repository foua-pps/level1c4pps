#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Martin.Raspaud

# Author(s):

#   Nina.Hakansson  a001865
#   Adam.Dybbroe
#   Martin.Raspaud

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

"""Script to make seviri level1c in PPS-format with pytroll"""

import argparse
from level1c4pps.avhrr_gac import process_one_file


# -----------------------------------------------------------------------------
# Main:
if __name__ == "__main__":
    """ Create PPS-format level1c data
    From a list of hirt files hrit create a level1c file for pps.
    """
    parser = argparse.ArgumentParser(
        description=('Script to produce a PPS-level1c file for a '
                     'GAC.'))
    parser.add_argument('file', type=str,
                        help='GAC file to process')
    parser.add_argument('-o', '--out_dir', type=str, nargs='?',
                        required=False, default='.',
                        help="Output directory where to store level1c file.")

    options = parser.parse_args()
    process_one_file(options.file, options.out_dir)
