#!/usr/bin/env python3

import io
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from surface_exchange import (
    SurfaceExchange,
    climate_vectors_in_body_frame,
    climate_controls,
    interpolate_spherical_surface,
    read_exchange,
    surface_to_climate_exchange,
    write_aspect_structured,
    write_exchange,
)


class SurfaceExchangeTests(unittest.TestCase):
    def test_round_trip_and_climate_controls(self):
        longitude_axis = np.array([-135.0, -45.0, 45.0, 135.0])
        latitude_axis = np.array([-45.0, 45.0])
        longitude, latitude = np.meshgrid(longitude_axis, latitude_axis)
        precipitation = np.linspace(1.0e-6, 8.0e-6, longitude.size)
        exchange = SurfaceExchange(
            1.0,
            longitude.reshape(-1),
            latitude.reshape(-1),
            {
                "precipitation_rate": precipitation,
                "surface_temperature": np.full(longitude.size, 288.15),
                "surface_elevation": np.zeros(longitude.size),
                "ice_thickness": np.linspace(0.0, 1000.0, longitude.size),
                "basal_ice_velocity": np.linspace(0.0, 100.0, longitude.size),
            },
            {
                "precipitation_rate": "kg m-2 s-1",
                "surface_temperature": "K",
                "surface_elevation": "m",
                "ice_thickness": "m",
                "basal_ice_velocity": "m yr-1",
            },
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "exchange.cxe"
            write_exchange(path, exchange)
            restored = read_exchange(path)
        np.testing.assert_allclose(restored.longitude, exchange.longitude)
        np.testing.assert_allclose(restored.fields["precipitation_rate"], precipitation)
        erosion_strength, runoff = climate_controls(restored)
        self.assertTrue(np.all(np.isfinite(erosion_strength)))
        self.assertTrue(np.all(runoff >= 0.0))

    def test_big_endian_exchange_and_aspect_spherical_coordinates(self):
        longitude_axis = np.array([-135.0, -45.0, 45.0, 135.0])
        latitude_axis = np.array([-45.0, 45.0])
        longitude, latitude = np.meshgrid(longitude_axis, latitude_axis)
        exchange = SurfaceExchange(
            2.0,
            longitude.reshape(-1),
            latitude.reshape(-1),
            {"elevation_change": np.arange(longitude.size, dtype=float)},
            {"elevation_change": "m"},
            ">",
        )
        with tempfile.TemporaryDirectory() as directory:
            binary_path = Path(directory) / "exchange.cxe"
            text_path = Path(directory) / "aspect.txt"
            write_exchange(binary_path, exchange)
            restored = read_exchange(binary_path)
            write_aspect_structured(
                text_path, restored, restored.fields["elevation_change"]
            )
            text = text_path.read_text(encoding="utf-8")
        self.assertEqual(restored.byte_order, ">")
        self.assertIn("# POINTS: 6 4", text)
        coordinates = np.loadtxt(io.StringIO("\n".join(text.splitlines()[2:])))
        self.assertGreaterEqual(coordinates[:, 0].min(), 0.0)
        self.assertLessEqual(coordinates[:, 0].max(), 2.0 * np.pi)
        self.assertGreaterEqual(coordinates[:, 1].min(), 0.0)
        self.assertLessEqual(coordinates[:, 1].max(), np.pi)

    def test_climate_coordinates_follow_spin_axis(self):
        longitude = np.array([0.0, 90.0, 0.0])
        latitude = np.array([0.0, 0.0, 90.0])
        geographic = climate_vectors_in_body_frame(
            longitude, latitude, np.array([0.0, 0.0, 1.0])
        )
        np.testing.assert_allclose(
            geographic,
            np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            atol=1.0e-14,
        )

        pole_on_body_x = climate_vectors_in_body_frame(
            longitude, latitude, np.array([1.0, 0.0, 0.0])
        )
        np.testing.assert_allclose(
            pole_on_body_x[-1], [1.0, 0.0, 0.0], atol=1.0e-14
        )
        np.testing.assert_allclose(
            pole_on_body_x[0], [0.0, 1.0, 0.0], atol=1.0e-14
        )

    def test_weighted_spherical_interpolation_is_exact_at_source_points(self):
        source = np.eye(3)
        values = np.array([2.0, 4.0, 8.0])
        interpolated = interpolate_spherical_surface(source, values, source, 3)
        np.testing.assert_allclose(interpolated, values, atol=1.0e-14)

    def test_command_writes_climate_topography(self):
        longitude, latitude = np.meshgrid(
            np.array([-135.0, -45.0, 45.0, 135.0]),
            np.array([-45.0, 45.0]),
        )
        exchange = SurfaceExchange(
            1.0,
            longitude.reshape(-1),
            latitude.reshape(-1),
            {
                "precipitation_rate": np.full(longitude.size, 1.0e-6),
                "surface_temperature": np.full(longitude.size, 288.15),
                "surface_elevation": np.arange(longitude.size, dtype=float),
                "ice_thickness": np.zeros(longitude.size),
                "basal_ice_velocity": np.zeros(longitude.size),
            },
            {},
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            exchange_path = directory / "climate.cxe"
            topography_path = directory / "topography.txt"
            erosion_path = directory / "erosion.txt"
            runoff_path = directory / "runoff.txt"
            write_exchange(exchange_path, exchange)
            subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).with_name("surface_exchange.py")),
                    "climate-to-aspect",
                    str(exchange_path),
                    "--surface-topography",
                    str(topography_path),
                    "--erosion-strength",
                    str(erosion_path),
                    "--surface-runoff",
                    str(runoff_path),
                ],
                check=True,
            )
            values = np.loadtxt(topography_path, skiprows=2)[:, 2]
        self.assertEqual(values.size, 24)
        self.assertGreater(values.max(), values.min())

    def test_surface_return_can_extract_one_window_increment(self):
        longitude_axis = np.array([-135.0, -45.0, 45.0, 135.0])
        latitude_axis = np.array([-45.0, 45.0])
        longitude, latitude = np.meshgrid(longitude_axis, latitude_axis)
        exchange = SurfaceExchange(
            1.0,
            longitude.reshape(-1),
            latitude.reshape(-1),
            {"surface_elevation": np.zeros(longitude.size)},
            {},
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            climate_path = directory / "climate.cxe"
            reference_path = directory / "reference.csv"
            current_path = directory / "current.csv"
            output_path = directory / "increment.cxe"
            write_exchange(climate_path, exchange)
            header = "longitude_deg,latitude_deg,elevation_change_m\n"
            reference_lines = [header]
            current_lines = [header]
            for lon, lat in zip(exchange.longitude, exchange.latitude):
                reference_lines.append(f"{lon},{lat},2\n")
                current_lines.append(f"{lon},{lat},5\n")
            reference_path.write_text("".join(reference_lines), encoding="utf-8")
            current_path.write_text("".join(current_lines), encoding="utf-8")
            surface_to_climate_exchange(
                current_path,
                climate_path,
                output_path,
                reference_surface_csv=reference_path,
            )
            returned = read_exchange(output_path)
        np.testing.assert_allclose(returned.fields["elevation_change"], 3.0)


if __name__ == "__main__":
    unittest.main()
