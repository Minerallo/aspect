#!/usr/bin/env python3
"""Show CLIMBER-X ice, its Fastscape remap, and polar-wander feedback."""

import argparse
from pathlib import Path
import sys
import tempfile

import matplotlib.pyplot as plt
import numpy as np


adapter_directory = Path(__file__).resolve().parents[2] / "contrib" / "coupling" / "climberx_fastscape"
sys.path.insert(0, str(adapter_directory))
from surface_exchange import (  # noqa: E402
    read_exchange,
    spin_axis_from_history,
    surface_to_climate_exchange,
)


def read_named_table(filename):
    return np.genfromtxt(filename, delimiter=",", names=True)


parser = argparse.ArgumentParser()
parser.add_argument("climate_exchange", type=Path)
parser.add_argument("fastscape_output", type=Path)
parser.add_argument("--output", type=Path, default=Path("climberx_ice_feedback.png"))
arguments = parser.parse_args()

climate = read_exchange(arguments.climate_exchange)
surface_directory = arguments.fastscape_output / "fastscape_surface_evolution"
surface = read_named_table(surface_directory / "surface-00001.csv")
pole = read_named_table(arguments.fastscape_output / "true_polar_wander.csv")

with tempfile.TemporaryDirectory() as temporary_directory:
    geographic_file = Path(temporary_directory) / "geographic.cxe"
    spin_frame_file = Path(temporary_directory) / "spin-frame.cxe"
    surface_file = surface_directory / "surface-00002.csv"
    surface_to_climate_exchange(
        surface_file, arguments.climate_exchange, geographic_file
    )
    surface_to_climate_exchange(
        surface_file,
        arguments.climate_exchange,
        spin_frame_file,
        spin_axis_from_history(arguments.fastscape_output / "true_polar_wander.csv"),
    )
    geographic_return = read_exchange(geographic_file)
    spin_frame_return = read_exchange(spin_frame_file)
return_difference = (
    spin_frame_return.fields["elevation_change"]
    - geographic_return.fields["elevation_change"]
)

figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
axes = axes.reshape(-1)
common = dict(cmap="Blues", vmin=0)

source_plot = axes[0].scatter(
    climate.longitude % 360,
    climate.latitude,
    c=climate.fields["ice_thickness"],
    s=15,
    **common,
)
axes[0].set(title="CLIMBER-X native exchange", xlabel="longitude (degrees)", ylabel="latitude (degrees)")
figure.colorbar(source_plot, ax=axes[0], label="grounded ice thickness (m)")

remapped_plot = axes[1].scatter(
    surface["longitude_deg"] % 360,
    surface["latitude_deg"],
    c=surface["ice_thickness_m"],
    s=13,
    **common,
)
axes[1].set(title="Field sampled by Fastscape", xlabel="longitude (degrees)", ylabel="latitude (degrees)")
figure.colorbar(remapped_plot, ax=axes[1], label="ice thickness (m)")

axes[2].plot(
    pole["time_years"],
    90.0-pole["pole_latitude_degrees"],
    "o-",
    label="spin-axis displacement",
)
axes[2].set(
    title="Ice-load feedback to polar wander",
    xlabel="time (years)",
    ylabel="angular distance from initial pole (degrees)",
)
axes[2].grid(alpha=0.25)
load_ratio = pole["effective_ice_load_kg_m2"] / pole["rigid_ice_load_kg_m2"]
load_axis = axes[2].twinx()
load_axis.plot(pole["time_years"], load_ratio, "s--", color="tab:purple", label="effective / rigid ice load")
load_axis.set_ylabel("self-gravity-corrected load ratio", color="tab:purple")
load_axis.tick_params(axis="y", colors="tab:purple")

return_limit = max(np.max(np.abs(return_difference)), 1.0e-12)
return_plot = axes[3].scatter(
    climate.longitude % 360,
    climate.latitude,
    c=return_difference,
    s=15,
    cmap="coolwarm",
    vmin=-return_limit,
    vmax=return_limit,
)
axes[3].set(
    title="Polar-wander change returned to CLIMBER-X",
    xlabel="longitude (degrees)",
    ylabel="latitude (degrees)",
)
figure.colorbar(return_plot, ax=axes[3], label="topography increment difference (m)")

for axis in (axes[0], axes[1], axes[3]):
    axis.set(xlim=(0, 360), ylim=(-90, 90))
    axis.grid(alpha=0.15)

figure.suptitle("CLIMBER-X ice distribution enters Fastscape erosion and true polar wander")
figure.savefig(arguments.output, dpi=180)
