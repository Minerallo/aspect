#!/usr/bin/env python3
"""Generate smooth Norway-like valley-glacier forcing fields."""

from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
LENGTH = 160_000.0
WIDTH = 100_000.0
NX = 17
NY = 11


def glacier_weight(x: float, y: float) -> float:
    """Return a trunk glacier plus two headward tributaries."""
    trunk = math.exp(-((y - 50_000.0) / 11_000.0) ** 2)
    tributary_offset = 0.28 * max(80_000.0 - x, 0.0)
    north = 0.72 * math.exp(
        -((y - (50_000.0 + tributary_offset)) / 8_000.0) ** 2
    )
    south = 0.72 * math.exp(
        -((y - (50_000.0 - tributary_offset)) / 8_000.0) ** 2
    )
    longitudinal = max(math.sin(math.pi * x / LENGTH), 0.0) ** 0.45
    return longitudinal * max(trunk, north, south)


def write_field(filename: str, column: str, value) -> None:
    lines = [f"# POINTS: {NX} {NY}", f"x y {column}"]
    for iy in range(NY):
        y = WIDTH * iy / (NY - 1)
        for ix in range(NX):
            x = LENGTH * ix / (NX - 1)
            lines.append(f"{x:.0f} {y:.0f} {value(x, y):.8g}")
    (ROOT / filename).write_text("\n".join(lines) + "\n")


write_field(
    "norway_ice_thickness.txt",
    "ice_thickness",
    lambda x, y: 1_400.0 * glacier_weight(x, y),
)
write_field(
    "norway_basal_ice_velocity.txt",
    "basal_ice_velocity",
    lambda x, y: ((35.0 + 145.0 * x / LENGTH) / 3.0) * glacier_weight(x, y),
)
write_field("zero_ice_thickness.txt", "ice_thickness", lambda _x, _y: 0.0)
write_field("zero_ice_velocity.txt", "basal_ice_velocity", lambda _x, _y: 0.0)
