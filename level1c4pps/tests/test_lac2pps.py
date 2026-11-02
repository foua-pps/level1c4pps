#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Unit tests for the lac2pps_lib module."""

import datetime as dt
import inspect
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import netCDF4
import numpy as np
import xarray as xr
from satpy import Scene

import level1c4pps.lac2pps_lib as lac2pps
from level1c4pps import save_data


TEST_DIR = os.path.dirname(__file__)
TIROSN_LEVEL1B_FILE = os.path.join(TEST_DIR, "NSS.GHRR.TN.D80003.S1147.E1332.B0630506.GC")
METOPB_LEVEL1B_FILE = os.path.join(TEST_DIR, "ESR.LHRR.M1.D16087.S2023.E2037.B01828628.BN")
SCANLINE_TIMES = np.array(["2009-07-01T12:16", "2009-07-01T12:27"], dtype="datetime64[ms]")
SATPY_ANGLE_NAMES = ("solar_zenith_angle", "sensor_zenith_angle", "solar_azimuth_angle",
                     "sensor_azimuth_angle", "sun_sensor_azimuth_difference_angle")


def _make_swath_array(name, values, dims=("y", "x"), **attrs):
    """Make an array along the scanlines, shaped like the ones the avhrr_l1b_gaclac reader loads."""
    return xr.DataArray(values, dims=dims, coords={"acq_time": ("y", SCANLINE_TIMES)}, attrs={"name": name, **attrs})


def _make_channel(name, values, wavelength, platform_name="noaa19"):
    """Make a channel carrying what the avhrr_l1b_gaclac reader gives every channel."""
    return _make_swath_array(name, values,
                             wavelength=wavelength,
                             start_time=dt.datetime(2009, 7, 1, 12, 16),
                             end_time=dt.datetime(2009, 7, 1, 12, 27),
                             platform_name=platform_name,
                             sensor="avhrr-3")


def _make_scene(platform_name="noaa19"):
    """Make a scene shaped like one the avhrr_l1b_gaclac reader loads."""
    scene = Scene()
    scene["4"] = _make_channel("4", [[280.0, 281.0], [282.0, 283.0]], [10.3, 10.8, 11.3, "um"], platform_name)
    scene["1"] = _make_channel("1", [[10.0, 11.0], [12.0, 13.0]], [0.58, 0.63, 0.68, "um"], platform_name)
    scene["latitude"] = _make_swath_array("latitude", [[60.0, 61.0], [62.0, 63.0]])
    scene["longitude"] = _make_swath_array("longitude", [[15.0, 16.0], [17.0, 18.0]])
    for name in SATPY_ANGLE_NAMES:
        scene[name] = _make_swath_array(name, [[10.0, 20.0], [30.0, 40.0]])
    scene["qual_flags"] = _make_swath_array("qual_flags", [[1, 0, 0, 0, 0, 0, 0], [2, 0, 0, 0, 0, 0, 0]],
                                            dims=("y", "num_flags"))
    return scene


class TestProcessScene(unittest.TestCase):
    """Test converting an already loaded scene to PPS level1c."""

    @staticmethod
    def _write(scene, **kwargs):
        """Convert the scene with the writer stubbed; return the file name and the writer's arguments by name."""
        with mock.patch("level1c4pps.lac2pps_lib.save_data") as writer:
            filename = lac2pps.process_scene(scene, out_path="/out", **kwargs)
        return filename, inspect.signature(save_data).bind(*writer.call_args.args, **writer.call_args.kwargs).arguments

    def _written_scene(self, scene):
        """Convert the scene with the writer stubbed and return the scene that reached the writer."""
        _, written = self._write(scene)
        return written["scene"]

    def _pps_identity(self, scene, channel):
        """Return the image name and tag the channel carries when it reaches the writer."""
        attrs = self._written_scene(scene)[channel].attrs
        return attrs["name"], attrs.get("id_tag")

    def _pps_identity_of_added_channel(self, channel, values, wavelength):
        """Add the channel to the shared scene and return the PPS identity it reaches the writer with."""
        scene = _make_scene()
        scene[channel] = _make_channel(channel, values, wavelength)
        return self._pps_identity(scene, channel)

    def test_the_scene_is_written_under_the_pps_file_name(self):
        """The file name is what PPS uses to find a pass, so it must follow the PPS pattern."""
        expected = os.path.join("/out", "S_NWC_avhrr_noaa19_12345_20090701T1216000Z_20090701T1227000Z.nc")
        filename, written = self._write(_make_scene(), orbit_n=12345)
        self.assertEqual((filename, written["filename"]), (expected, expected))

    def test_channel_1_reaches_the_writer_as_the_06_micron_reflectance_image(self):
        """PPS finds a channel by its image name and tag, never by the AVHRR channel name."""
        self.assertEqual(self._pps_identity(_make_scene(), "1"), ("image1", "ch_r06"))

    def test_channel_3_of_older_platforms_reaches_the_writer_as_the_37_micron_image(self):
        """AVHRR/1 and AVHRR/2 call their 3.7 micron channel 3; PPS must still find it as image5."""
        identity = self._pps_identity_of_added_channel("3", [[290.0, 291.0], [292.0, 293.0]], [3.55, 3.74, 3.93, "um"])
        self.assertEqual(identity, ("image5", "ch_tb37"))

    def test_channel_2_reaches_the_writer_as_the_09_micron_reflectance_image(self):
        """A LAC or FRAC file must look to PPS exactly like a GAC one from the same instrument."""
        identity = self._pps_identity_of_added_channel("2", [[20.0, 21.0], [22.0, 23.0]], [0.725, 0.8625, 1.0, "um"])
        self.assertEqual(identity, ("image2", "ch_r09"))

    def test_channel_3a_reaches_the_writer_as_the_16_micron_reflectance_image(self):
        """AVHRR/3 measures 1.6 micron in channel 3a by day; PPS knows that channel as image6."""
        identity = self._pps_identity_of_added_channel("3a", [[5.0, 6.0], [7.0, 8.0]], [1.58, 1.61, 1.64, "um"])
        self.assertEqual(identity, ("image6", "ch_r16"))

    def test_channel_3b_reaches_the_writer_as_the_37_micron_image(self):
        """AVHRR/3 measures 3.7 micron in channel 3b; PPS knows it as image5, as for channel 3 on older platforms."""
        identity = self._pps_identity_of_added_channel("3b", [[250.0, 251.0], [252.0, 253.0]], [3.55, 3.74, 3.93, "um"])
        self.assertEqual(identity, ("image5", "ch_tb37"))

    def test_channel_5_reaches_the_writer_as_the_12_micron_image(self):
        """Channel 5 is the 12 micron split-window channel; PPS knows it as image4."""
        identity = self._pps_identity_of_added_channel("5", [[270.0, 271.0], [272.0, 273.0]], [11.5, 12.0, 12.5, "um"])
        self.assertEqual(identity, ("image4", "ch_tb12"))

    def test_channel_5_of_an_avhrr_1_platform_does_not_reach_the_writer(self):
        """AVHRR/1 has four channels, yet pygac gives it a channel 5 that is a copy of channel 4.

        PPS would take that channel for a real 12 micron channel.
        """
        for platform_name in ("tirosn", "noaa6", "noaa8", "noaa10"):
            with self.subTest(platform_name=platform_name):
                scene = _make_scene(platform_name)
                scene["5"] = _make_channel("5", [[280.0, 281.0], [282.0, 283.0]], [11.5, 12.0, 12.5, "um"],
                                           platform_name)
                self.assertNotIn("5", self._written_scene(scene))

    def test_a_noaa_10_pass_without_channel_5_still_reaches_the_writer(self):
        """Once the reader stops offering NOAA-10 a channel 5, its passes must still convert instead of crashing."""
        self.assertIn("4", self._written_scene(_make_scene("noaa10")))

    def test_channel_4_reaches_the_writer_as_the_11_micron_image(self):
        """Channel 4 is the 11 micron window channel every AVHRR carries; PPS knows it as image3."""
        self.assertEqual(self._pps_identity(_make_scene(), "4"), ("image3", "ch_tb11"))

    def test_latitude_and_longitude_reach_the_writer_as_lat_and_lon(self):
        """PPS reads its geolocation from variables called lat and lon, and from nothing else."""
        written_scene = self._written_scene(_make_scene())
        self.assertEqual([name in written_scene for name in ("lat", "lon", "latitude", "longitude")],
                         [True, True, False, False])

    def test_the_angles_reach_the_writer_under_the_names_pps_reads_them_by(self):
        """PPS knows the sun and satellite geometry only under its own five angle names."""
        written_scene = self._written_scene(_make_scene())
        pps_names = ("sunzenith", "satzenith", "sunazimuth", "satazimuth", "azimuthdiff")
        self.assertEqual(([name in written_scene for name in pps_names],
                          [name in written_scene for name in SATPY_ANGLE_NAMES]),
                         ([True] * 5, [False] * 5))

    def test_the_sun_zenith_angle_reaches_the_writer_with_its_pps_tag_and_cf_standard_name(self):
        """PPS selects an angle by its tag; the standard name tells any other reader what it holds."""
        written_scene = self._written_scene(_make_scene())
        attrs = written_scene["sunzenith"].attrs
        self.assertEqual((attrs.get("id_tag"), attrs.get("standard_name")), ("sunzenith", "solar_zenith_angle"))

    def test_every_scanline_reaches_the_writer_with_its_acquisition_time(self):
        """PPS times each line separately; without the timestamps a pass has only a start and an end."""
        written_scene = self._written_scene(_make_scene())
        np.testing.assert_array_equal(written_scene["scanline_timestamps"].values, SCANLINE_TIMES)

    def test_the_scanline_times_reach_the_writer_only_as_scanline_timestamps(self):
        """PPS reads the time of each scanline from scanline_timestamps, as in GAC files.

        The reader also hangs that time on every dataset as an acq_time coordinate.
        Written as it is, the file carries the scanline times twice: once more as an
        acq_time variable that every image names among its coordinates.
        """
        written_scene = self._written_scene(_make_scene())
        carrying_acq_time = sorted(data.attrs["name"] for data in written_scene.values() if "acq_time" in data.coords)
        self.assertEqual(carrying_acq_time, [])

    def test_the_quality_flags_reach_the_writer_with_the_pps_tag_and_name_of_the_gac_flags(self):
        """PPS reads the pygac quality flags under the same tag and name whether a pass is GAC or LAC."""
        attrs = self._written_scene(_make_scene())["qual_flags"].attrs
        self.assertEqual((attrs.get("id_tag"), attrs.get("long_name")), ("qual_flags", "pygac quality flags"))

    def test_the_header_names_lac2pps_as_its_source(self):
        """The source attribute is how a PPS user tells which converter made the file."""
        _, written = self._write(_make_scene())
        self.assertEqual((written["header_attrs"] or {}).get("source"), "lac2pps.py")

    def test_the_writer_uses_the_netcdf_engine_it_is_asked_for(self):
        """pygac-fdr writes with the engine its configuration names, so that choice must reach the writer."""
        _, written = self._write(_make_scene(), engine="netcdf4")
        self.assertEqual(written["engine"], "netcdf4")

    def test_the_header_carries_the_orbit_number(self):
        """PPS keeps the orbit number in the file header as well as in the file name."""
        _, written = self._write(_make_scene(), orbit_n=12345)
        self.assertEqual((written["header_attrs"] or {}).get("orbit_number"), 12345)

    def test_the_header_carries_what_every_level1c4pps_avhrr_header_carries(self):
        """PPS reads platform, sensor, times and provenance from a GAC file's header; a LAC file must offer the same.

        Conventions is not listed: the CF writer adds it to the file itself.
        """
        _, written = self._write(_make_scene(), orbit_n=12345)
        gac_header_keys = {"date_created", "end_time", "history", "instrument", "orbit_number", "platform",
                           "platform_name", "sensor", "source", "start_time", "version_level1c4pps",
                           "version_level1c4pps_satpy"}
        self.assertEqual(gac_header_keys - set(written["header_attrs"] or {}), set())


class TestProcessOneFile(unittest.TestCase):
    """Test converting a level 1b file to a PPS level1c file on disk."""

    def test_a_four_channel_level1b_file_becomes_a_pps_file_with_its_four_images(self):
        """The whole chain, from the reader to a file on disk, on the oldest AVHRR there is: TIROS-N."""
        with tempfile.TemporaryDirectory() as out_dir:
            filename = lac2pps.process_one_file(TIROSN_LEVEL1B_FILE, out_path=out_dir,
                                                reader_kwargs={"tle_dir": TEST_DIR, "tle_name": "TLE_tirosn.txt"})
            with netCDF4.Dataset(filename) as pps_file:
                images = sorted(name for name in pps_file.variables if name.startswith("image"))
        self.assertEqual(images, ["image1", "image2", "image3", "image5"])

    def test_an_avhrr_3_lac_file_becomes_a_pps_file_with_its_six_images(self):
        """AVHRR/3 splits channel 3 into 3a and 3b and adds channel 5; PPS needs all six images from a LAC pass."""
        with tempfile.TemporaryDirectory() as out_dir:
            filename = lac2pps.process_one_file(METOPB_LEVEL1B_FILE, out_path=out_dir,
                                                reader_kwargs={"tle_dir": TEST_DIR, "tle_name": "TLE_metopb.txt"})
            with netCDF4.Dataset(filename) as pps_file:
                images = sorted(name for name in pps_file.variables if name.startswith("image"))
        self.assertEqual(images, ["image1", "image2", "image3", "image4", "image5", "image6"])


class TestLac2ppsScript(unittest.TestCase):
    """Test the lac2pps command-line script."""

    def test_the_script_writes_one_pps_file_for_a_level1b_file(self):
        """Operators run the converters from the command line, so the conversion must be reachable from there."""
        script = os.path.join(TEST_DIR, os.pardir, os.pardir, "bin", "lac2pps.py")
        with tempfile.TemporaryDirectory() as out_dir:
            run = subprocess.run([sys.executable, script, TIROSN_LEVEL1B_FILE, "-o", out_dir, "-td", TEST_DIR,
                                  "-tn", "TLE_tirosn.txt"], capture_output=True)
            written = [name for name in os.listdir(out_dir) if name.startswith("S_NWC_avhrr_")]
        self.assertEqual((run.returncode, len(written)), (0, 1))
