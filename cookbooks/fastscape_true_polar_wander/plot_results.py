#!/usr/bin/env python3
"""Plot the pole path and the climate forcing displaced over the surface."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_named_table(filename):
    return np.genfromtxt(filename, delimiter=",", names=True)


parser = argparse.ArgumentParser()
parser.add_argument("output_directory", type=Path)
parser.add_argument("--output", type=Path, default=Path("true_polar_wander_feedback.png"))
arguments = parser.parse_args()

pole = read_named_table(arguments.output_directory / "true_polar_wander.csv")
surface_directory = arguments.output_directory / "fastscape_surface_evolution"
surface_1 = read_named_table(surface_directory / "surface-00001.csv")
surface_2 = read_named_table(surface_directory / "surface-00002.csv")

figure, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
axes[0].plot(pole["pole_longitude_degrees"], pole["pole_latitude_degrees"], "o-", label="evolving pole")
axes[0].plot(pole["equilibrium_longitude_degrees"], pole["equilibrium_latitude_degrees"], "x--", label="maximum inertia axis")
for time, longitude, latitude in zip(pole["time_years"], pole["pole_longitude_degrees"], pole["pole_latitude_degrees"]):
    axes[0].annotate(f"{time:g} yr", (longitude, latitude), xytext=(4, 4), textcoords="offset points", fontsize=8)
axes[0].set(xlim=(0, 360), ylim=(-90, 90), xlabel="body-fixed longitude (degrees)", ylabel="latitude (degrees)", title="Spin-axis evolution")
axes[0].grid(alpha=0.25)
axes[0].legend(fontsize=8)

common = dict(s=12, cmap="Blues", vmin=0)
first = axes[1].scatter(surface_1["longitude_deg"] % 360, surface_1["latitude_deg"], c=surface_1["ice_thickness_m"], **common)
axes[1].set(title="Ice sampled after 1,000 years", xlabel="longitude (degrees)", ylabel="latitude (degrees)", xlim=(0, 360), ylim=(-90, 90))
figure.colorbar(first, ax=axes[1], label="ice thickness (m)")

difference = surface_2["surface_runoff_factor"] - surface_1["surface_runoff_factor"]
limit = max(np.max(np.abs(difference)), 1e-12)
second = axes[2].scatter(surface_2["longitude_deg"] % 360, surface_2["latitude_deg"], c=difference, s=12, cmap="coolwarm", vmin=-limit, vmax=limit)
axes[2].set(title="Runoff change from moving climate coordinates", xlabel="longitude (degrees)", ylabel="latitude (degrees)", xlim=(0, 360), ylim=(-90, 90))
figure.colorbar(second, ax=axes[2], label="runoff factor change")

for axis in axes[1:]:
    axis.grid(alpha=0.15)

figure.suptitle("ASPECT–FastScape–CLIMBER-X true-polar-wander feedback test")
figure.savefig(arguments.output, dpi=180)
