#!/usr/bin/env python3
"""Check incision, sediment routing, and unloading in the Norway-like run."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def read(output: str, step: int) -> dict[str, np.ndarray]:
    path = ROOT / output / "fastscape_surface_evolution" / f"surface-{step:05d}.csv"
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {
        field: np.asarray([float(row[field]) for row in rows])
        for field in rows[0]
        if field not in {"bedrock_lithology"}
    }


loaded = read("output-fastscape-norway-glacial", 20)
unloaded = read("output-fastscape-norway-glacial", 22)
control = read("output-fastscape-norway-control", 20)

assert np.allclose(loaded["surface_x_m"], control["surface_x_m"])
assert np.allclose(loaded["surface_y_m"], control["surface_y_m"])

ice = loaded["ice_thickness_m"] >= 100.0
fjord = ice & (np.abs(loaded["surface_y_m"] - 50_000.0) <= 12_000.0)
incision = control["elevation_m"] - loaded["elevation_m"]
mean_fjord_incision = float(np.mean(incision[fjord]))
maximum_last_step_glacial_erosion = float(np.max(loaded["glacial_erosion_m"]))
maximum_sediment_flux = float(np.max(loaded["sediment_flux_m3_per_year"]))
maximum_marine_deposit = float(np.max(loaded["sediment_thickness_m"]))
maximum_rebound_velocity = float(
    np.max(unloaded["regional_ice_load_velocity_m_per_year"])
)

assert mean_fjord_incision > 10.0
assert maximum_last_step_glacial_erosion > 0.05
assert maximum_sediment_flux > 0.0
assert maximum_rebound_velocity > 0.01
assert np.max(unloaded["ice_thickness_m"]) == 0.0

print(f"mean fjord incision relative to control: {mean_fjord_incision:.2f} m")
print(f"maximum glacial erosion in last step:    {maximum_last_step_glacial_erosion:.3f} m")
print(f"maximum sediment flux:                   {maximum_sediment_flux:.2f} m^3/yr")
print(f"maximum marine deposit thickness:        {maximum_marine_deposit:.3f} m")
print(f"maximum unloading rebound velocity:      {maximum_rebound_velocity:.4f} m/yr")
