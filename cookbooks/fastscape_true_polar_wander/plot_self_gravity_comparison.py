#!/usr/bin/env python3
"""Compare rigid and degree-two self-gravitating ice-load responses."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_pole_history(output_directory):
    return np.genfromtxt(
        output_directory / "true_polar_wander.csv",
        delimiter=",",
        names=True,
    )


parser = argparse.ArgumentParser()
parser.add_argument("rigid_output", type=Path)
parser.add_argument("self_gravity_output", type=Path)
parser.add_argument("--output", type=Path, default=Path("self_gravity_comparison.png"))
arguments = parser.parse_args()

rigid = read_pole_history(arguments.rigid_output)
self_gravity = read_pole_history(arguments.self_gravity_output)

figure, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
axes[0].plot(rigid["time_years"], 90-rigid["pole_latitude_degrees"], "o-", label="rigid ice load")
axes[0].plot(self_gravity["time_years"], 90-self_gravity["pole_latitude_degrees"], "o-", label="degree-two load response")
axes[0].set(xlabel="time (years)", ylabel="angular distance from initial pole (degrees)", title="Polar motion")
axes[0].grid(alpha=0.25)
axes[0].legend()

response_ratio = self_gravity["effective_ice_load_kg_m2"] / self_gravity["rigid_ice_load_kg_m2"]
axes[1].plot(self_gravity["time_years"], response_ratio, "o-", color="tab:purple")
axes[1].axhline(0.1, linestyle="--", color="0.4", label="1 + fluid load Love number")
axes[1].set(xlabel="time (years)", ylabel="effective / rigid load-tensor norm", title="Self-gravitational compensation", ylim=(0, 1.05))
axes[1].grid(alpha=0.25)
axes[1].legend()

figure.suptitle("Degree-two self-gravity reduces ice-load-driven polar wander")
figure.savefig(arguments.output, dpi=180)
