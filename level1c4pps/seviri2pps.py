#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Nina.Hakansson, Adam.Dybbroe, Martin.Raspaud

# Author(s):

#   Nina.Hakansson  
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

# This program was developed by CMSAF to be used for the processing of
# CLAAS3. 

"""Script to make seviri level1c in PPS-format with pytroll"""

import os
import sys
import numpy as np
import xarray as xr
import dask.array as da
from glob import glob
import time
from datetime import datetime
from satpy.scene import Scene
from trollsift.parser import globify, Parser
from pyorbital.astronomy import get_alt_az, sun_zenith_angle
from pyorbital.orbital import get_observer_look
from calibration_coefs import get_calibration_for_time, CALIB_MODE
import pyresample
class UnexpectedSatPyVersion(Exception):
    pass


BANDNAMES = ['VIS006', 'VIS008', 'IR_016', 'IR_039',
             'IR_087', 'IR_108', 'IR_120', 
             'IR_134', 'IR_097', 'WV_062', 'WV_073']
PPS_TAGNAMES = {'VIS006': 'ch_r06', 
                'VIS008': 'ch_r09', 
                'IR_016': 'ch_r16', 
                'IR_039': 'ch_tb37',
                'IR_087': 'ch_tb85', 
                'IR_108': 'ch_tb11', 
                'IR_120': 'ch_tb12', 
                'IR_134': 'ch_tb133', 
                'IR_097': 'ch_tb97', 
                'WV_062': 'ch_tb67', 
                'WV_073': 'ch_tb73'}
PLATFORM_SHORTNAMES = {"MSG1": "Meteosat-8",
                       "MSG2": "Meteosat-9",
                       "MSG3": "Meteosat-10",
                       "MSG4": "Meteosat-11"}

# H-000-MSG3__-MSG3________-IR_120___-000003___-201410051115-__:
hrit_file_pattern = '{rate:1s}-000-{hrit_format:_<6s}-{platform_shortname:_<12s}-{channel:_<8s}_-{segment:_<9s}-{start_time:%Y%m%d%H%M}-__'
p__ = Parser(hrit_file_pattern)

def make_azidiff_angle(sata, suna, fill):
    """ Calculate azimuth difference angle """
    #np.ma.mod(np.ma.abs(sunaz - sataz), 180) same as?
    daz = sata - suna
    daz[daz < 0] = -1 * daz[daz < 0]
    daz[daz > 360] = daz[daz > 360] - 360
    daz[daz > 180] = 360 - daz[daz > 180]
    # fix nodata
    #daz[suna == fill] = fill
    #daz[sata == fill] = fill
    return daz

def process_one_scan(tslot_files, out_path, 
                     process_buggy_satellite_zenith_angles=False):
    """ Make level 1c files in PPS-format """
    tic = time.time()
    image_num = 0 # name of first dataset is image0
    #if len(tslot_files) != 8 * len(BANDNAMES) + 2:
    #    raise Exception("Some data is missing")
    platform_shortname = p__.parse(
        os.path.basename(tslot_files[0]))['platform_shortname']
    start_time = p__.parse(
        os.path.basename(tslot_files[0]))['start_time']
    platform_name = PLATFORM_SHORTNAMES[platform_shortname]
    #Load channel data for one scene and set some attributes
    coefs = get_calibration_for_time(platform=platform_shortname, 
                                     time=start_time)

    scn_ = Scene(
        reader='seviri_l1b_hrit',
        filenames=tslot_files,
        reader_kwargs={
            'calib_mode': CALIB_MODE,
            'ext_calib_coefs': coefs})
    scn_.attrs['platform_name'] = platform_name

    #SEVIRI data only
    if scn_.attrs['sensor'] == {'seviri'}: 
        sensor = 'seviri'
        scn_.load(BANDNAMES)        
    for band in BANDNAMES:
        idtag = PPS_TAGNAMES[band]
        scn_[band].attrs['id_tag'] = idtag
        scn_[band].attrs['description'] = 'SEVIRI ' + str(band)
        scn_[band].attrs['sun_earth_distance_correction_applied'] = 'False'
        scn_[band].attrs['sun_earth_distance_correction_factor'] = 1.0    
        scn_[band].attrs['sun_zenith_angle_correction_applied'] = 'False'
        scn_[band].attrs['name'] = "image{:d}".format(image_num)
        scn_[band].attrs['coordinates'] = 'lon lat'
        image_num += 1   

    #correct area  
    area_corr = pyresample.geometry.AreaDefinition(
        'seviri-corrected',
        'Corrected SEVIRI L1.5 grid (since Dec 2017)',
        'geosmsg',
        {'a': 6378169.00, 'b': 6356583.80, 'h': 35785831.0, 
         'lon_0': 0.0, 'proj': 'geos', 'units': 'm'},
        3712, 3712,
        (5567248.28340708, 5570248.686685662, 
         -5570248.686685662, -5567248.28340708)
    )
    if not scn_['IR_108'].attrs['georef_offset_corrected']:
        scn_ = scn_.resample(area_corr)
        print(scn_['IR_108'].attrs['georef_offset_corrected'])
    
    #import pdb;pdb.set_trace()
    #Set som header attributes:
    scn_.attrs['platform'] = platform_name
    scn_.attrs['instrument'] = sensor.upper()
    scn_.attrs['source'] = "seviri2pps.py"
    scn_.attrs['orbit_number'] = "99999"
    #scn_.attrs['orbit'] = "99999"
    nowutc = datetime.utcnow()
    scn_.attrs['date_created'] = nowutc.strftime("%Y-%m-%dT%H:%M:%SZ")        
    #Find lat/lon data 
    irch = scn_['IR_108']
    lons, lats = irch.attrs['area'].get_lonlats()
    lons[lons>360] = -999.0
    lons[lons<-360] = -999.0
    lats[lats>360] = -999.0
    lats[lats<-360] = -999.0

    #Find angles data
    sunalt, suna = get_alt_az(
        irch.attrs['start_time'], *irch.attrs['area'].get_lonlats())
    suna = np.rad2deg(suna)
    sunz = sun_zenith_angle(
        irch.attrs['start_time'], *irch.attrs['area'].get_lonlats())

    # if:
    #   Buggy data is requested buggy data is prepared!
    # elif:
    #   1) get_observer_look() gives wrong answer ...
    #   ... for satellite altitude in m. AND
    #   2) get_observer_look() gives correct answer ...
    #   ....  for satellite altitude in km. AND
    #   3) Satellite altitude is m.:
    #    => Satellite alltitude need to be converted to km.
    # else:
    #    => There have been updates to SatPy and this script 
    #       need to be modified.
    if process_buggy_satellite_zenith_angles:
        print(" Making buggy satellite zenith angels on purpose!")
        sata, satel = get_observer_look(
            irch.attrs['navigation']['satellite_actual_longitude'],
            irch.attrs['navigation']['satellite_actual_latitude'],
            irch.attrs['navigation']['satellite_actual_altitude'],
            irch.attrs['start_time'],
            lons, lats, 0)  
    elif (get_observer_look(0, 0, 36000*1000, 
                          datetime.utcnow(), np.array([16]), 
                          np.array([58]), np.array([0]))[1]>30 and
        get_observer_look(0, 0, 36000, 
                          datetime.utcnow(), np.array([16]), 
                          np.array([58]), np.array([0]))[1]<23 and 
        irch.attrs['satellite_altitude']>38000):
        sata, satel = get_observer_look(
            irch.attrs['navigation']['satellite_actual_longitude'],
            irch.attrs['navigation']['satellite_actual_latitude'],
            0.001*irch.attrs['navigation']['satellite_actual_altitude'],
            irch.attrs['start_time'],
            lons, lats, 0)
    else:
        raise UnexpectedSatpyVersion(
            "You might have a newer version of satpy/pyorbital that"
            "handles units. In that case the m => km conversion might"
            "be unneeded and wrong.")
            
    satz = 90 - satel
    azidiff = make_azidiff_angle(sata, suna, -32767)
    #Add lat/lon  and angles datasets to the scen object
    my_coords = scn_['IR_108'].coords
    my_coords['time'] = irch.attrs['start_time']
    scn_['lat'] = xr.DataArray(
        da.from_array(lats, chunks=(53, 3712)), 
        dims=['y','x'], 
        coords={'y': scn_['IR_108']['y'], 'x': scn_['IR_108']['x']})
    scn_['lat'].attrs['long_name'] = 'latitude coordinate'
    scn_['lat'].attrs['standard_name'] = 'latitude'
    scn_['lat'].attrs['units'] = 'degrees_north'
    scn_['lat'].attrs['start_time'] = irch.attrs['start_time']
    scn_['lat'].attrs['end_time'] = irch.attrs['end_time']
    scn_['lon'] = xr.DataArray(
        da.from_array(lons, chunks=(53, 3712)),  
        dims=['y','x'], 
        coords={'y': scn_['IR_108']['y'], 'x': scn_['IR_108']['x']})
    scn_['lon'].attrs['long_name'] = 'longitude coordinate'
    scn_['lon'].attrs['standard_name'] = 'longitude'
    scn_['lon'].attrs['units'] = 'degrees_east'
    scn_['lon'].attrs['start_time'] = irch.attrs['start_time']
    scn_['lon'].attrs['end_time'] = irch.attrs['end_time']
    #sunzenith
    scn_['sunzenith'] = xr.DataArray(
        da.from_array(sunz[:,:], chunks=(53, 3712)),  
        dims=['y','x'], coords=my_coords)
    scn_['sunzenith'].attrs['id_tag'] = 'sunzenith'
    scn_['sunzenith'].attrs['long_name'] = 'sun zenith angle'
    scn_['sunzenith'].attrs['standard_name'] = 'solar_zenith_angle'
    scn_['sunzenith'].attrs['valid_range'] = [0, 18000]
    scn_['sunzenith'].attrs['name'] = "image{:d}".format(image_num)
    image_num +=1
    #satzenith
    scn_['satzenith'] = xr.DataArray(
        da.from_array(satz[:,:], chunks=(53, 3712)),            
        dims=['y','x'], coords=my_coords)
    scn_['satzenith'].attrs['id_tag'] = 'satzenith'
    scn_['satzenith'].attrs['long_name'] = 'satellite zenith angle'
    scn_['satzenith'].attrs['standard_name'] = 'platform_zenith_angle'
    scn_['satzenith'].attrs['valid_range'] = [0, 9000]
    scn_['satzenith'].attrs['name'] = "image{:d}".format(image_num)
    image_num +=1
    #azidiff
    scn_['azimuthdiff'] = xr.DataArray(
        da.from_array(azidiff[:,:], chunks=(53, 3712)),
        dims=['y','x'], coords=my_coords)
    scn_['azimuthdiff'].attrs['id_tag'] = 'azimuthdiff'
    scn_['azimuthdiff'].attrs['standard_name'] = (
        'angle_of_rotation_from_solar_azimuth_to_platform_azimuth')
    scn_['azimuthdiff'].attrs['long_name'] = 'azimuth difference angle'
    scn_['azimuthdiff'].attrs['valid_range'] = [0, 18000]
    scn_['azimuthdiff'].attrs['name'] = "image{:d}".format(image_num)
    image_num +=1
    for angle in ['azimuthdiff', 'satzenith', 'sunzenith']:
        scn_[angle].attrs['units'] =  'degree'
        for attr in irch.attrs.keys():
            if attr in ["start_time",
                        "end_time",
                        "navigation",
                        "georef_offset_corrected",
                        "projection"
            ]:
                scn_[angle].attrs[attr] = irch.attrs[attr]

    #Get filename
    start_time = scn_['IR_108'].attrs['start_time']
    end_time = scn_['IR_108'].attrs['end_time']
    filename = os.path.join(
        out_path, 
        "S_NWC_seviri_{:s}_{:s}_{:s}Z_{:s}Z.nc".format(
            platform_name.lower().replace('-',''),
            "99999",
            start_time.strftime('%Y%m%dT%H%M%S%f')[:-5],
            end_time.strftime('%Y%m%dT%H%M%S%f')[:-5]))
        
    #Encoding for channels
    save_info = {}
    for band in BANDNAMES:
        idtag = PPS_TAGNAMES[band]
        name = scn_[band].attrs['name']
        scn_[band].attrs.pop('area',None)
        # Add time coordinate. To make cfwriter aware that we want 3D data.
        my_coords = scn_[band].coords
        my_coords['time'] = irch.attrs['start_time']

        if 'tb' in idtag:
            save_info[name] = {'dtype': 'int16', 
                               'scale_factor':0.01, 
                               '_FillValue': -32767, 
                               'zlib': True,
                               'complevel': 4,
                               'add_offset': 273.15 }
        else:
            save_info[name] = {'dtype': 'int16', 
                               'scale_factor':0.01, 
                               'zlib': True,
                               'complevel': 4,
                               '_FillValue': -32767, 
                               'add_offset': 0.0 }
    #Encoding for angles and lat/lon       
    for name in ['image11', 'image12', 'image13']:    
        save_info[name] = {
            'dtype': 'int16', 
            'scale_factor':0.01, 
            'zlib': True,
            'complevel': 4,
            '_FillValue': -32767, 
            'add_offset': 0.0 }

    for name in ['lon', 'lat']:
        save_info[name] = {'dtype': 'float32',    'zlib': True,
                           'complevel': 4, '_FillValue': -999.0}
    header_attrs = scn_.attrs.copy()    
    header_attrs['start_time'] = time.strftime(
        "%Y-%m-%d %H:%M:%S",
        irch.attrs['start_time'].timetuple())
    header_attrs['end_time'] = time.strftime(
        "%Y-%m-%d %H:%M:%S",
        irch.attrs['end_time'].timetuple())
    header_attrs['sensor'] = sensor.lower()
    header_attrs.pop('platform_name', None)
            

    scn_.save_datasets(writer='cf', 
                       filename=filename,  
                       header_attrs=header_attrs, 
                       engine='netcdf4', 
                       encoding=save_info,
                       include_lonlats=False,
                       pretty=True,
                       flatten_attrs=True,
                       exclude_attrs=['raw_metadata'])
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic)) #About 40 seconds 

def process_all_scans_in_dname(dname, out_path, ok_dates=None):
    """ Make level 1c files for all files in directory dname """
    fl_ = glob(os.path.join(dname, globify(hrit_file_pattern)))
    dates = [p__.parse(os.path.basename(p))['start_time'] for p in fl_]
    unique_dates = np.unique(dates).tolist()

    for uqdate in unique_dates:
        date_formated = uqdate.strftime("%Y%m%d%H%M")
        if ok_dates is not None and date_formated not in ok_dates.keys():
            print("Skipping date {date}".format(date=date_formated))
            continue
        # Every hour only:
        #if uqdate.minute != 0:
        #    continue
        tslot_files = [f for f in fl_ if p__.parse(
            os.path.basename(f))['start_time'] == uqdate]
        try:
            process_one_scan(tslot_files, out_path)
        except:
            pass


# -----------------------------------------------------------------------------
# Main:
if __name__ == "__main__":
    """ Create PPS-format level1c data 
    From a list of hirt files hrit create a level1c file for pps.
    """
    #python3 seviri2pps.py file1 file2 ... fileN -o output
    import argparse
    parser = argparse.ArgumentParser(
        description = ('Script to produce a PPS-level1c file for a list of '
                       'SEVIRI hrit files.'))
    parser.add_argument('files', metavar='fileN', type=str, nargs='+',
                        help='List of hrit files to process for one scan')
    parser.add_argument('-o', '--out_dir', type=str, nargs='?', 
                        required=False, 
                        help="Output directory where to store level1c file.")
    parser.add_argument('-b', '--buggy_satz', const=True, nargs='?', 
                        required=False,
                        help="Create buggy satellite zenith angle data")
    options = parser.parse_args()
    process_buggy_satellite_zenith_angles=False
    if options.buggy_satz:
        process_buggy_satellite_zenith_angles=True
    process_one_scan(options.files, options.out_dir, 
                     process_buggy_satellite_zenith_angles)



