#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path

import numpy as np

from surface_exchange import SurfaceExchange, climate_controls, read_exchange, write_exchange


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
            },
            {
                "precipitation_rate": "kg m-2 s-1",
                "surface_temperature": "K",
                "surface_elevation": "m",
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


if __name__ == "__main__":
    unittest.main()
