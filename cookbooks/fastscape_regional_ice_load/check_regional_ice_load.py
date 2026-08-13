#!/usr/bin/env python3
"""Check regional ice loading, unloading, and the disabled reference run."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def read_surface(output: str, step: int) -> list[dict[str, str]]:
    path = ROOT / output / "fastscape_surface_evolution" / f"surface-{step:05d}.csv"
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def mean(rows: list[dict[str, str]], field: str) -> float:
    return sum(float(row[field]) for row in rows) / len(rows)


def central_cells(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        row for row in rows
        if 25_000 <= float(row["surface_x_m"]) <= 75_000
        and 25_000 <= float(row["surface_y_m"]) <= 75_000
    ]


loaded = central_cells(read_surface("output-fastscape-regional-ice-load", 5))
unloaded_path = (
    ROOT
    / "output-fastscape-regional-ice-load"
    / "fastscape_surface_evolution"
    / "surface-00010.csv"
)
with unloaded_path.open(newline="") as stream:
    unloaded = list(csv.DictReader(stream))
unloaded_center = central_cells(unloaded)
reference = central_cells(
    read_surface("output-fastscape-regional-ice-load-no-response", 10)
)

mean_ice_thickness = mean(loaded, "ice_thickness_m")
equilibrium = -917.0 * mean_ice_thickness / 3300.0
expected_loaded = equilibrium * (1.0 - math.exp(-1.0))
expected_unloaded = expected_loaded * math.exp(-1.0)
loaded_displacement = mean(loaded, "regional_ice_load_displacement_m")
loaded_elevation = mean(loaded, "elevation_m")
unloaded_displacement = mean(
    unloaded_center, "regional_ice_load_displacement_m"
)
unloaded_elevation = mean(unloaded_center, "elevation_m")
unloaded_velocity = max(
    float(row["regional_ice_load_velocity_m_per_year"]) for row in unloaded
)
reference_displacement = mean(reference, "regional_ice_load_displacement_m")

assert math.isclose(loaded_displacement, expected_loaded, rel_tol=1e-10)
assert math.isclose(unloaded_displacement, expected_unloaded, rel_tol=1e-10)
assert math.isclose(loaded_elevation, loaded_displacement, abs_tol=1e-8)
assert math.isclose(unloaded_elevation, unloaded_displacement, abs_tol=1e-8)
assert unloaded_velocity > 0.0
assert abs(reference_displacement) < 1e-12

print(f"loaded displacement:   {loaded_displacement:.3f} m")
print(f"unloaded displacement: {unloaded_displacement:.3f} m")
print(f"rebound velocity:       {unloaded_velocity:.6f} m/yr")
print("disabled reference:     0 m")
