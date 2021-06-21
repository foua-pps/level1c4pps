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
#   Stephan Finkensieper <stephan.finkensieper@dwd.de>

# This program was developed by CMSAF to be used for the processing of
# CLAAS3.

"""Tools to convert SEVIRI hrit to PPS level-1c format."""


import os
import numpy as np
import xarray as xr
import dask.array as da
from glob import glob
import time
from datetime import datetime, timedelta
from satpy.scene import Scene
import satpy.utils
from satpy.readers.utils import remove_earthsun_distance_correction
from trollsift.parser import globify, Parser
from pyorbital.astronomy import get_alt_az, sun_zenith_angle
from pyorbital.orbital import get_observer_look

from level1c4pps.calibration_coefs import get_calibration_for_time, CALIB_MODE
from level1c4pps import make_azidiff_angle, get_encoding, compose_filename, update_angle_attributes


try:
    FileNotFoundError
except NameError:
    # Python 2
    FileNotFoundError = IOError


class UnexpectedSatpyVersion(Exception):
    """Exception if unexpected satpy version."""

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


HRIT_FILE_PATTERN = ('{rate:1s}-000-{hrit_format:_<6s}-'
                     '{platform_shortname:_<12s}-{channel:_<8s}_-'
                     '{segment:_<9s}-{start_time:%Y%m%d%H%M}-__')
# H-000-MSG3__-MSG3________-IR_120___-000003___-201410051115-__

NATIVE_FILE_PATTERN = ('{platform_shortname:4s}-{instr:4s}-'
                       'MSG{product_level:2d}-{base_algorithm_version:4s}-NA-'
                       '{end_time:%Y%m%d%H%M%S.%f}000Z')
# MSG4-SEVI-MSG15-0100-NA-20190409121243.927000000Z-20190409121300-1329370.nat
# MSG4-SEVI-MSG15-0201-NA-20190409121243.927000000Z.nat
# MSG4-SEVI-MSG15-1234-NA-20190409121243.927000000Z


def load_and_calibrate(filenames, apply_sun_earth_distance_correction, rotate):
    """Load and calibrate data.

    Uses inter-calibration coefficients from Meirink et al.

    Args:
        filenames: List of data files
        apply_sun_earth_distance_correction: If True, apply sun-earth-distance-
            correction to visible channels.
        rotate: Rotate image so that pixel (0, 0) is N-W.

    Returns:
        Satpy scene holding calibrated channels
    """
    # Parse filenames
    parser = SEVIRIFilenameParser()
    file_format, info = parser.parse(os.path.basename(filenames[0]))

    calib_coefs = get_calibration_for_time(
        platform=info['platform_shortname'],
        time=info['start_time']
    )
    scn_ = _create_scene(file_format, filenames, calib_coefs)
    _check_is_seviri_data(scn_)
    _load_bands(scn_, rotate)
    _update_scene_attrs(scn_, {'image_rotated': rotate})

    if not apply_sun_earth_distance_correction:
        remove_sun_earth_distance_correction(scn_)

    return scn_


def _create_scene(file_format, filenames, calib_coefs):
    return Scene(reader=file_format,
                 filenames=filenames,
                 reader_kwargs={'calib_mode': CALIB_MODE,
                                'ext_calib_coefs': calib_coefs})


def _check_is_seviri_data(scene):
    if not scene.attrs['sensor'] == {'seviri'}:
        raise ValueError('Not SEVIRI data')


def _load_bands(scene, rotate):
    scene.load(
        BANDNAMES,
        upper_right_corner=_get_upper_right_corner(rotate)
    )


def _get_upper_right_corner(rotation_flag):
    return 'NE' if rotation_flag else 'native'


def _update_scene_attrs(scene, attrs):
    scene.attrs.update(attrs)


def remove_sun_earth_distance_correction(scene):
    """Remove sun-earth-distance correction from visible channels.

    This is required because a downstream CLAAS module (CPP) does not recognize
    the "sun_earth_distance_correction_applied" attribute and applies the
    correction anyway.
    """
    for band in BANDNAMES:
        is_vis = scene[band].attrs['calibration'] == 'reflectance'
        correction_applied = scene[band].attrs.get(
            'sun_earth_distance_correction_applied', False
        )
        if is_vis and correction_applied:
            scene[band] = remove_earthsun_distance_correction(scene[band])


def get_lonlats(dataset):
    """Get lat/lon coordinates."""
    lons, lats = dataset.attrs['area'].get_lonlats()
    lons[np.fabs(lons) > 360] = np.nan
    lats[np.fabs(lats) > 90] = np.nan
    return lons, lats


def get_solar_angles(scene, lons, lats):
    """Compute solar angles.

    Compute angles for each scanline using their acquisition time to account for
    the earth's rotation over the course of one scan.

    Returns:
        Solar azimuth angle, Solar zenith angle in degrees

    """
    suna = np.full(lons.shape, np.nan)
    sunz = np.full(lons.shape, np.nan)
    mean_acq_time = get_mean_acq_time(scene)
    for line, acq_time in enumerate(mean_acq_time.values):
        if np.isnat(acq_time):
            continue
        _, suna_line = get_alt_az(acq_time, lons[line, :], lats[line, :])
        suna_line = np.rad2deg(suna_line)
        suna[line, :] = suna_line
        sunz[line, :] = sun_zenith_angle(acq_time, lons[line, :], lats[line, :])
    return suna, sunz


def get_satellite_angles(dataset, lons, lats):
    """Compute satellite angles.

    Returns:
        Satellite azimuth angle, Satellite zenith angle in degrees

    """
    sat_lon, sat_lat, sat_alt = satpy.utils.get_satpos(dataset)

    # Double check that pyorbital/satpy behave as expected (satpy returning
    # altitude in meters and pyorbital expecting km).
    #
    # if:
    #   1) get_observer_look() gives wrong answer ...
    #   ... for satellite altitude in m. AND
    #   2) get_observer_look() gives correct answer ...
    #   ....  for satellite altitude in km. AND
    #   3) Satellite altitude is m.:
    #    => Satellite alltitude need to be converted to km.
    # else:
    #    => There have been updates to SatPy and this script
    #       need to be modified.
    if not (get_observer_look(0, 0, 36000*1000,
                              datetime.utcnow(), np.array([16]),
                              np.array([58]), np.array([0]))[1] > 30 and
            get_observer_look(0, 0, 36000,
                              datetime.utcnow(), np.array([16]),
                              np.array([58]), np.array([0]))[1] < 23 and
            sat_alt > 38000):
        raise UnexpectedSatpyVersion(
            'Unexpected handling of satellite altitude in pyorbital/'
            'satpy. Conversion to km is probably unneeded and wrong.')

    # Convert altitude from meters to kilometers, as expected by the
    # current version of pyorbital
    sat_alt *= 0.001

    # Compute angles
    sata, satel = get_observer_look(
        sat_lon,
        sat_lat,
        sat_alt,
        dataset.attrs['start_time'],
        lons, lats, 0)
    satz = 90 - satel

    return sata, satz


def set_attrs(scene):
    """Set global and band attributes."""
    # Global
    scene.attrs['platform'] = scene['IR_108'].attrs['platform_name']
    scene.attrs['instrument'] = 'SEVIRI'
    scene.attrs['source'] = "seviri2pps.py"
    scene.attrs['orbit_number'] = 99999
    nowutc = datetime.utcnow()
    scene.attrs['date_created'] = nowutc.strftime("%Y-%m-%dT%H:%M:%SZ")

    # For each band
    for image_num, band in enumerate(BANDNAMES):
        idtag = PPS_TAGNAMES[band]
        scene[band].attrs['id_tag'] = idtag
        scene[band].attrs['description'] = 'SEVIRI ' + str(band)
        scene[band].attrs['wavelength'] = [scene[band].attrs['wavelength'].min,
                                           scene[band].attrs['wavelength'].central,
                                           scene[band].attrs['wavelength'].max]
        if 'sun_earth_distance_correction_factor' not in scene[band].attrs:
            scene[band].attrs['sun_earth_distance_correction_applied'] = False
            scene[band].attrs['sun_earth_distance_correction_factor'] = 1.0
        scene[band].attrs['sun_zenith_angle_correction_applied'] = False
        scene[band].attrs['name'] = "image{:d}".format(image_num)

        # Cosmetics
        for attr in ['orbital_parameters', 'satellite_longitude',
                     'satellite_latitude', 'satellite_altitude',
                     'platform_name', 'sensor', 'georef_offset_corrected']:
            scene[band].attrs.pop(attr, None)


def get_mean_acq_time(scene):
    """Compute mean scanline acquisition time over all bands."""
    dtype = scene['IR_108'].coords['acq_time'].dtype

    # Convert timestamps to float to facilitate averaging. Caveat: NaT is
    # not converted to NaN, but to -9.22E18. So we have to set these elements
    # to NaN manually
    acq_times = []
    for band in BANDNAMES:
        acq_time = scene[band].coords['acq_time'].drop_vars(['acq_time'])
        is_nat = np.isnat(acq_time.values)
        acq_time = acq_time.astype(int).where(np.logical_not(is_nat))
        acq_times.append(acq_time)

    # Compute average over all bands (skip NaNs)
    acq_times = xr.concat(acq_times, 'bands')
    return acq_times.mean(dim='bands', skipna=True).astype(dtype)


def update_coords(scene):
    """Update band coordinates."""
    mean_acq_time = get_mean_acq_time(scene)
    for band in BANDNAMES:
        # Override channel-specific scanline timestamps with mean acquisition
        # time. The differences are not very large and the resulting nc file
        # is much simpler.
        scene[band]['acq_time'] = mean_acq_time

        # Remove area, set lat/lon as coordinates
        scene[band].attrs['coordinates'] = 'lon lat'
        area = scene[band].attrs.pop('area', None)
        if area:
            scene.attrs['area'] = area

        # Add time coordinate to make cfwriter aware that we want 3D data
        scene[band].coords['time'] = scene[band].attrs['start_time']


def add_ancillary_datasets(scene, lons, lats, sunz, satz, azidiff,
                           chunks=(512, 3712)):
    """Add ancillary datasets to the scene.

    Args:
        lons: Longitude coordinates
        lats: Latitude coordinates
        sunz: Solar zenith angle
        satz: Satellite zenith angle
        azidiff: Absolute azimuth difference angle
        chunks: Chunksize

    """
    start_time = scene['IR_108'].attrs['start_time']
    end_time = scene['IR_108'].attrs['end_time']
    angle_coords = scene['IR_108'].coords

    # Latitude
    scene['lat'] = xr.DataArray(
        da.from_array(lats, chunks=chunks),
        dims=['y', 'x'],
        coords={'y': scene['IR_108']['y'], 'x': scene['IR_108']['x']})
    scene['lat'].attrs['long_name'] = 'latitude coordinate'
    scene['lat'].attrs['standard_name'] = 'latitude'
    scene['lat'].attrs['units'] = 'degrees_north'
    scene['lat'].attrs['start_time'] = start_time
    scene['lat'].attrs['end_time'] = end_time

    # Longitude
    scene['lon'] = xr.DataArray(
        da.from_array(lons, chunks=chunks),
        dims=['y', 'x'],
        coords={'y': scene['IR_108']['y'], 'x': scene['IR_108']['x']})
    scene['lon'].attrs['long_name'] = 'longitude coordinate'
    scene['lon'].attrs['standard_name'] = 'longitude'
    scene['lon'].attrs['units'] = 'degrees_east'
    scene['lon'].attrs['start_time'] = start_time
    scene['lon'].attrs['end_time'] = end_time

    # Sunzenith
    scene['sunzenith'] = xr.DataArray(
        da.from_array(sunz[:, :], chunks=chunks),
        dims=['y', 'x'], coords=angle_coords)

    # Satzenith
    scene['satzenith'] = xr.DataArray(
        da.from_array(satz[:, :], chunks=chunks),
        dims=['y', 'x'], coords=angle_coords)

    # Azidiff
    scene['azimuthdiff'] = xr.DataArray(
        da.from_array(azidiff[:, :], chunks=chunks),
        dims=['y', 'x'], coords=angle_coords)

    # Update the attributes
    update_angle_attributes(scene, band=scene['IR_108'])


def add_proj_satpos(scene):
    """Add information on projection and satellite position."""
    orb = scene['IR_108'].attrs['orbital_parameters']

    # Area extent
    try:
        # Traditionally a, b was included in proj_dict
        param_a = scene.attrs['area'].proj_dict['a']
        param_b = scene.attrs['area'].proj_dict['b']
    except KeyError:
        # But sometimes b is missing:
        param_a = scene.attrs['area'].crs.ellipsoid.semi_major_metre
        param_b = scene.attrs['area'].crs.ellipsoid.semi_minor_metre

    scene.attrs.update({
        'projection': 'geos',
        'projection_semi_major_axis': param_a,
        'projection_semi_minor_axis': param_b,
        'projection_longitude': orb['projection_longitude'],
        'projection_latitude': orb['projection_latitude'],
        'projection_altitude': orb['projection_altitude']
    })
    scene['projection_area_extent'] = xr.DataArray(
        np.array([scene.attrs['area'].area_extent]),
        dims=('time', 'corners'),
        attrs={'units': 'meters',
               'long_name': 'Projection area extent',
               'comment': 'Area extent (lower left x, lower left y, '
                          'upper right x, upper right y) in projection '
                          'coordinates'}
    )
    scene['georef_offset_corrected'] = xr.DataArray(
        np.array([scene['IR_108'].attrs['georef_offset_corrected']],
                 dtype='int8'),
        dims=('time', ),
        attrs={'long_name': 'Georeferencing offset correction flag',
               'flag_values': np.array([0, 1], dtype='int8'),
               'flag_meanings': 'georef_offset_present georef_offset_corrected',
               'comment': 'Until December 2017, SEVIRI L1.5 data is shifted '
                          'by 1.5km SSP North and West against the nominal '
                          'GEOS projection. Since December 2017 this offset '
                          'has been corrected by EUMETSAT. This flag indicates '
                          'if the correction has been applied. If not, the area '
                          'extent is adjusted so that the projection always '
                          'matches the data.'}
    )

    # Satellite position
    scene['satellite_longitude'] = xr.DataArray(
        [orb['satellite_actual_longitude']],
        dims=('time',),
        attrs={'units': 'degree_east',
               'long_name': 'Actual satellite longitude',
               'comment': 'Not to be confused with projection longitude'}
    )
    scene['satellite_latitude'] = xr.DataArray(
        [orb['satellite_actual_latitude']],
        dims=('time',),
        attrs={'units': 'degree_north',
               'long_name': 'Actual satellite latitude',
               'comment': 'Not to be confused with projection latitude'}
    )
    scene['satellite_altitude'] = xr.DataArray(
        [orb['satellite_actual_altitude']],
        dims=('time',),
        attrs={'units': 'meters',
               'long_name': 'Actual satellite altitude',
               'comment': 'Relative and perpendicular to the surface of the '
                          'earth ellipsoid. Not to be confused with projection '
                          'altitude.'}
    )


def get_encoding_seviri(scene):
    """Get netcdf encoding for all datasets."""
    # Bands
    chunks = (1, 512, 3712)
    encoding = get_encoding(scene,
                            bandnames=BANDNAMES,
                            pps_tagnames=PPS_TAGNAMES,
                            chunks=chunks)

    # Time
    acq_units = scene.attrs['start_time'].strftime(
        'milliseconds since %Y-%m-%d %H:%M')
    encoding['acq_time'] = {'units': acq_units,
                            'calendar': 'standard',
                            '_FillValue': -9999.0}
    encoding['time'] = {'units': 'days since 2004-01-01 00:00',
                        'calendar': 'standard',
                        '_FillValue': None,
                        'chunksizes': [1]}

    return encoding


def get_header_attrs(scene):
    """Get global netcdf attributes."""
    header_attrs = scene.attrs.copy()
    for key in ['sensor', 'area']:
        header_attrs.pop(key, None)
    return header_attrs


def _drop_seconds_milliseconds(dtime):
    return dtime.replace(second=0, microsecond=0)


def set_nominal_scan_time(dataset):
    """Set time attributes in the dataset to the nominal scan time."""
    dataset = dataset.copy()
    start_time = _drop_seconds_milliseconds(
        dataset.attrs['start_time'])
    end_time = start_time + timedelta(minutes=15)
    dataset.attrs['start_time'] = start_time
    dataset.attrs['end_time'] = end_time
    return dataset


class SEVIRIFilenameParser:
    """SEVIRI filename parser."""

    default_formats = [
        {'name': 'seviri_l1b_hrit',
         'pattern': HRIT_FILE_PATTERN,
         'full_match': True},
        {'name': 'seviri_l1b_native',
         'pattern': NATIVE_FILE_PATTERN,
         'full_match': False},
    ]

    def __init__(self, formats=None):
        self.formats = formats or self.default_formats
        self.format_names = [fmt['name'] for fmt in self.formats]

    def parse(self, filename):
        """Parse the given filename.

        Returns:
            File format, filename info
        """
        for fmt in self.formats:
            parser = Parser(fmt['pattern'])
            try:
                info = parser.parse(filename, full_match=fmt['full_match'])
                break
            except ValueError:
                continue
        else:
            raise ValueError(
                'Filename {} doesn\'t match any of the supported '
                'formats {}.'.format(filename, self.format_names)

            )
        return fmt['name'], self._postproc(fmt['name'], info)

    def _postproc(self, format_name, info):
        """Postprocess filename info."""
        if format_name == 'seviri_l1b_native':
            return self._postproc_native(info)
        elif format_name == 'seviri_l1b_hrit':
            return self._postproc_hrit(info)
        else:
            raise NotImplementedError

    def _postproc_native(self, info):
        """Postprocess Native filename info."""
        # Only end time available, derive start time.
        end_time = info['end_time']
        quarter = end_time.minute // 15
        info['start_time'] = end_time.replace(
            minute=quarter * 15,
            second=0,
            microsecond=0
        )
        return info

    def _postproc_hrit(self, parsed):
        """Postprocess HRIT filename info."""
        return parsed


def process_one_scan(tslot_files, out_path, rotate=True, engine='h5netcdf',
                     use_nominal_time_in_filename=False,
                     apply_sun_earth_distance_correction=True):
    """Make level 1c files in PPS-format."""
    for fname in tslot_files:
        if not os.path.isfile(fname):
            raise FileNotFoundError('No such file: {}'.format(fname))

    tic = time.time()
    scn_ = load_and_calibrate(
        tslot_files,
        apply_sun_earth_distance_correction=apply_sun_earth_distance_correction,
        rotate=rotate
    )

    # Find lat/lon data
    lons, lats = get_lonlats(scn_['IR_108'])

    # Compute angles
    suna, sunz = get_solar_angles(scn_, lons=lons, lats=lats)
    sata, satz = get_satellite_angles(scn_['IR_108'], lons=lons, lats=lats)
    azidiff = make_azidiff_angle(sata, suna)

    # Update coordinates
    update_coords(scn_)

    # Add ancillary datasets to the scene
    add_ancillary_datasets(scn_, lons=lons, lats=lats, sunz=sunz, satz=satz,
                           azidiff=azidiff)
    add_proj_satpos(scn_)

    # Set attributes. This changes SEVIRI band names to PPS band names.
    set_attrs(scn_)

    # Write datasets to netcdf
    ir108_for_filename = scn_['IR_108']
    if use_nominal_time_in_filename:
        ir108_for_filename = set_nominal_scan_time(ir108_for_filename)
    filename = compose_filename(
        scene=scn_,
        out_path=out_path,
        instrument='seviri',
        band=ir108_for_filename
    )
    scn_.save_datasets(writer='cf',
                       filename=filename,
                       header_attrs=get_header_attrs(scn_),
                       engine=engine,
                       encoding=get_encoding_seviri(scn_),
                       unlimited_dims=['time'],
                       include_lonlats=False,
                       pretty=True,
                       flatten_attrs=True,
                       exclude_attrs=['raw_metadata'])
    print("Saved file {:s} after {:3.1f} seconds".format(
        os.path.basename(filename),
        time.time()-tic))  # About 40 seconds
    return filename


def process_all_scans_in_dname(dname, out_path, ok_dates=None, rotate=False):
    """Make level 1c files for all files in directory dname."""
    parser = Parser(HRIT_FILE_PATTERN)
    fl_ = glob(os.path.join(dname, globify(HRIT_FILE_PATTERN)))
    dates = [parser.parse(os.path.basename(p))['start_time'] for p in fl_]
    unique_dates = np.unique(dates).tolist()
    for uqdate in unique_dates:
        date_formated = uqdate.strftime("%Y%m%d%H%M")
        if ok_dates is not None and date_formated not in ok_dates.keys():
            print("Skipping date {date}".format(date=date_formated))
            continue
        # Every hour only:
        # if uqdate.minute != 0:
        #    continue
        tslot_files = [f for f in fl_ if parser.parse(
            os.path.basename(f))['start_time'] == uqdate]
        try:
            process_one_scan(tslot_files, out_path, rotate=rotate)
        except Exception:
            pass
