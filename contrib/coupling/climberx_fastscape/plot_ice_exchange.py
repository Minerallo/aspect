#!/usr/bin/env python3
"""Plot ice thickness and basal sliding speed from two native exchanges."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from surface_exchange import SurfaceExchange, read_exchange


def field_on_grid(exchange: SurfaceExchange, name: str) -> tuple[np.ndarray, ...]:
    longitude = np.unique(exchange.longitude)
    latitude = np.unique(exchange.latitude)
    field = np.full((latitude.size, longitude.size), np.nan)
    longitude_index = np.searchsorted(longitude, exchange.longitude)
    latitude_index = np.searchsorted(latitude, exchange.latitude)
    field[latitude_index, longitude_index] = exchange.fields[name]
    return longitude, latitude, field


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("first", type=Path)
    parser.add_argument("second", type=Path)
    parser.add_argument("--first-label", default="first ice model")
    parser.add_argument("--second-label", default="second ice model")
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()

    exchanges = [read_exchange(arguments.first), read_exchange(arguments.second)]
    labels = [arguments.first_label, arguments.second_label]
    thickness_maximum = max(
        np.nanmax(exchange.fields["ice_thickness"]) for exchange in exchanges
    )
    speed_maximum = max(
        np.nanmax(exchange.fields["basal_ice_velocity"]) for exchange in exchanges
    )

    figure, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    for row, (exchange, label) in enumerate(zip(exchanges, labels)):
        longitude, latitude, thickness = field_on_grid(exchange, "ice_thickness")
        _, _, speed = field_on_grid(exchange, "basal_ice_velocity")
        thickness_plot = axes[row, 0].pcolormesh(
            longitude, latitude, thickness, shading="auto", vmin=0, vmax=thickness_maximum
        )
        speed_plot = axes[row, 1].pcolormesh(
            longitude,
            latitude,
            np.log10(1.0 + speed),
            shading="auto",
            vmin=0,
            vmax=np.log10(1.0 + speed_maximum),
        )
        axes[row, 0].set_title(f"{label}: ice thickness")
        axes[row, 1].set_title(f"{label}: basal sliding speed")
        for axis in axes[row]:
            axis.set_xlabel("longitude (degrees)")
            axis.set_ylabel("latitude (degrees)")
        figure.colorbar(thickness_plot, ax=axes[row, 0], label="metres")
        figure.colorbar(
            speed_plot,
            ax=axes[row, 1],
            label="log10(1 + speed in metres per year)",
        )

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.suptitle("Dynamic ice fields passed to surface-process coupling")
    figure.savefig(arguments.output, dpi=180)


if __name__ == "__main__":
    main()
