#!/usr/bin/env python3

import io
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


if __name__ == "__main__":
    unittest.main()
