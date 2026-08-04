#!/usr/bin/env python3
"""Plot fields and changes from one closed physical feedback cycle."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

COUPLING = Path(__file__).resolve().parents[2] / "contrib/coupling/climberx_fastscape"
sys.path.insert(0, str(COUPLING))
from surface_exchange import read_exchange  # noqa: E402


def grid(exchange, values):
    longitude = np.unique(exchange.longitude)
    latitude = np.unique(exchange.latitude)
    return (
        longitude,
        latitude,
        np.asarray(values).reshape(latitude.size, longitude.size),
    )


def surface_column(path: Path, name: str):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return (
        np.array([float(row["longitude_deg"]) for row in rows]),
        np.array([float(row["latitude_deg"]) for row in rows]),
        np.array([float(row[name]) for row in rows]),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--feedback", required=True, type=Path)
    parser.add_argument("--topography", required=True, type=Path)
    parser.add_argument("--surface", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    baseline = read_exchange(arguments.baseline)
    feedback = read_exchange(arguments.feedback)
    returned = read_exchange(arguments.topography)
    panels = (
        (
            baseline.fields["precipitation_rate"] * 31536000.0,
            "Initial precipitation",
            "kg m⁻² yr⁻¹",
        ),
        (baseline.fields["ice_thickness"], "Dynamic ice thickness", "m"),
        (baseline.fields["basal_ice_velocity"], "Basal ice speed", "m yr⁻¹"),
        (
            returned.fields["elevation_change"],
            "Surface change returned to climate",
            "m",
        ),
        (
            feedback.fields["surface_temperature"]
            - baseline.fields["surface_temperature"],
            "Climate temperature response",
            "K",
        ),
        (
            feedback.fields["ice_thickness"] - baseline.fields["ice_thickness"],
            "Ice-thickness response",
            "m",
        ),
    )
    figure, axes = plt.subplots(3, 2, figsize=(13, 10), constrained_layout=True)
    for axis, (values, title, unit) in zip(axes.flat, panels):
        longitude, latitude, field = grid(baseline, values)
        color_map = "coolwarm" if np.min(values) < 0 else "viridis"
        image = axis.pcolormesh(
            longitude, latitude, field, shading="auto", cmap=color_map
        )
        axis.set(title=title, xlabel="longitude (degrees)", ylabel="latitude (degrees)")
        figure.colorbar(image, ax=axis, label=unit)
    longitude, latitude, glacial_erosion = surface_column(
        arguments.surface, "glacial_erosion_m"
    )
    active = glacial_erosion > 0
    axes[1, 0].scatter(
        longitude[active],
        latitude[active],
        s=4,
        c="white",
        alpha=0.5,
        label="glacial erosion",
    )
    if np.any(active):
        axes[1, 0].legend(loc="lower left")
    figure.suptitle(
        "Closed climate–ice–surface processes–solid Earth feedback cycle", fontsize=15
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output, dpi=180)


if __name__ == "__main__":
    main()
