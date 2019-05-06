#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Nina.Hakansson

# Author(s):

#   Nina.Hakansson  

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

import os
from glob import glob
from seviri2pps import process_all_scans_in_dname


if __name__ == "__main__":
    #For testing:
    my_tar_file = "local_disk/DATA_MISC/SEVIRI/20100208_15_MSG2.tar"
    my_temp_dir = "local_disk/temp/SEVIRI"
    out_path = "local_disk/DATA_MISC/SEVIRI/"
    sdate = "20100208"
    stime = "201002081715"

    tar_command = "tar -xvf %s -C %s"%(my_tar_file, my_temp_dir)
    files = glob("%s/%s/%s/*"%(my_temp_dir,sdate,stime)) 
    for filename in files:
        bzip2_command = "bzip2 -d %s"%(filename)
    dirl = glob("%s/%s/%s/"%(my_temp_dir,sdate,stime))  

    for dirname in dirl:
        print("Directory: ", dirname)
        if not os.path.isdir(dirname):
            print("Can't find directory.")
            continue
        try:
            process_all_scans_in_dname(dirname, out_path)
        except:
            raise

            
