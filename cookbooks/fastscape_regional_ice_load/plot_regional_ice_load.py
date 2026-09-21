#!/usr/bin/env python3
"""Plot regional ice-load subsidence, unloading rebound, and a control run."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def read_surface(output: str, step: int) -> dict[str, np.ndarray]:
    path = ROOT / output / "fastscape_surface_evolution" / f"surface-{step:05d}.csv"
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {
        name: np.array([float(row[name]) for row in rows])
        for name in (
            "surface_x_m",
            "surface_y_m",
            "regional_ice_load_displacement_m",
            "regional_ice_load_velocity_m_per_year",
        )
    }


def central_mean(surface: dict[str, np.ndarray], field: str) -> float:
    center = (
        (surface["surface_x_m"] >= 25_000)
        & (surface["surface_x_m"] <= 75_000)
        & (surface["surface_y_m"] >= 25_000)
        & (surface["surface_y_m"] <= 75_000)
    )
    return float(np.mean(surface[field][center]))


active_output = "output-fastscape-regional-ice-load"
control_output = "output-fastscape-regional-ice-load-no-response"
loaded = read_surface(active_output, 5)
unloaded = read_surface(active_output, 10)
times = np.arange(1, 11, dtype=float) * 1000.0
active_history = np.array(
    [central_mean(read_surface(active_output, step),
                  "regional_ice_load_displacement_m") for step in range(1, 11)]
)
control_history = np.array(
    [central_mean(read_surface(control_output, step),
                  "regional_ice_load_displacement_m") for step in range(1, 11)]
)

figure, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
limit = max(abs(loaded["regional_ice_load_displacement_m"]).max(),
            abs(unloaded["regional_ice_load_displacement_m"]).max())
for axis, surface, title in zip(
    axes[:2],
    (loaded, unloaded),
    ("After 5,000 years of loading", "After 5,000 years of rebound"),
):
    image = axis.scatter(
        surface["surface_x_m"] / 1000.0,
        surface["surface_y_m"] / 1000.0,
        c=surface["regional_ice_load_displacement_m"],
        s=12,
        cmap="coolwarm",
        vmin=-limit,
        vmax=limit,
    )
    axis.set_title(title)
    axis.set_xlabel("x (kilometres)")
    axis.set_ylabel("y (kilometres)")
    axis.set_aspect("equal")
figure.colorbar(image, ax=axes[:2], label="bedrock displacement (metres)")

axes[2].plot(times, active_history, "o-", label="loading then unloading")
axes[2].plot(times, control_history, "--", label="response disabled")
axes[2].axvline(5000, color="black", linestyle=":", label="ice removed")
axes[2].set_xlabel("time (years)")
axes[2].set_ylabel("mean displacement beneath ice block (metres)")
axes[2].set_title("Subsidence and rebound")
axes[2].legend()
axes[2].grid(alpha=0.25)

figure.suptitle("Regional ice-load response combined with Fastscape C++")
figure.savefig(ROOT / "regional_ice_load_response.png", dpi=180)
