#!/usr/bin/env python3
"""Render a cinematic Norway-like ice, erosion, and unloading animation."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".render-cache"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation, colors
import numpy as np
import pyvista as pv


OUTPUT = ROOT / "output-fastscape-norway-glacial" / "fastscape_surface_evolution"
VERTICAL_EXAGGERATION = 12.0
FRAMES_BETWEEN_STEPS = 3


def read_surface(step: int) -> dict[str, np.ndarray]:
    path = OUTPUT / f"surface-{step:05d}.csv"
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    fields = (
        "surface_x_m",
        "surface_y_m",
        "elevation_m",
        "ice_thickness_m",
        "basal_ice_velocity_m_per_year",
        "glacial_erosion_m",
        "sediment_thickness_m",
        "sediment_flux_m3_per_year",
        "regional_ice_load_velocity_m_per_year",
    )
    return {
        field: np.asarray([float(row[field]) for row in rows], dtype=float)
        for field in fields
    }


def regular_grid(data: dict[str, np.ndarray], stride: int) -> dict[str, np.ndarray]:
    """Put cell-centred output on a downsampled regular grid."""
    x_values = np.unique(data["surface_x_m"])
    y_values = np.unique(data["surface_y_m"])
    nx = len(x_values)
    ny = len(y_values)
    order = np.lexsort((data["surface_x_m"], data["surface_y_m"]))
    gridded = {}
    for field, values in data.items():
        gridded[field] = values[order].reshape(ny, nx)[::stride, ::stride]
    return gridded


def interpolate(first: dict[str, np.ndarray], second: dict[str, np.ndarray],
                fraction: float) -> dict[str, np.ndarray]:
    return {
        field: (1.0 - fraction) * first[field] + fraction * second[field]
        for field in first
    }


def pyvista_surface(frame: dict[str, np.ndarray]) -> pv.StructuredGrid:
    """Construct the terrain and attach every rendered diagnostic."""
    x = frame["surface_x_m"] / 1000.0
    y = frame["surface_y_m"] / 1000.0
    z = frame["elevation_m"] / 1000.0 * VERTICAL_EXAGGERATION
    grid = pv.StructuredGrid(x, y, z)
    for field, values in frame.items():
        grid.point_data[field] = values.ravel(order="F")
    return grid


def render(output: Path, poster: Path, stride: int, frames_per_second: int) -> None:
    surfaces = [regular_grid(read_surface(step), stride) for step in range(1, 23)]
    movie_frames: list[tuple[dict[str, np.ndarray], float]] = []
    for index in range(len(surfaces) - 1):
        for subframe in range(FRAMES_BETWEEN_STEPS):
            fraction = subframe / FRAMES_BETWEEN_STEPS
            movie_frames.append(
                (interpolate(surfaces[index], surfaces[index + 1], fraction),
                 5000.0 * (index + 1 + fraction))
            )
    movie_frames.append((surfaces[-1], 110000.0))

    all_elevation = np.concatenate([frame["elevation_m"].ravel() for frame in surfaces])
    elevation_norm = colors.Normalize(np.percentile(all_elevation, 1),
                                      np.percentile(all_elevation, 99))
    erosion_max = max(float(np.max(frame["glacial_erosion_m"])) for frame in surfaces)

    figure = plt.figure(figsize=(16, 9), facecolor="#08131e")
    axis = figure.add_axes((0.035, 0.07, 0.79, 0.88), projection="3d")
    axis.set_facecolor("#08131e")
    status_axis = figure.add_axes((0.82, 0.08, 0.16, 0.84), facecolor="#0d2030")
    status_axis.set_axis_off()

    def draw(frame_number: int) -> None:
        frame, time_years = movie_frames[frame_number]
        mesh = pyvista_surface(frame)
        points = mesh.points.reshape(frame["elevation_m"].shape + (3,), order="F")
        x = points[:, :, 0]
        y = points[:, :, 1]
        terrain_z = points[:, :, 2]
        ice_z = terrain_z + frame["ice_thickness_m"] / 1000.0 * VERTICAL_EXAGGERATION
        has_ice = frame["ice_thickness_m"] >= 100.0

        axis.clear()
        status_axis.clear()
        status_axis.set_axis_off()
        axis.set_facecolor("#08131e")

        terrain_colors = plt.colormaps["gist_earth"](
            elevation_norm(frame["elevation_m"])
        )
        terrain_colors[..., 3] = 1.0
        axis.plot_surface(x, y, terrain_z, facecolors=terrain_colors,
                          linewidth=0, antialiased=False, shade=True)

        if np.any(has_ice):
            ice_plot = np.where(has_ice, ice_z, np.nan)
            speed = frame["basal_ice_velocity_m_per_year"]
            speed_norm = colors.Normalize(0.0, max(float(np.max(speed)), 1.0))
            ice_colors = plt.colormaps["Blues_r"](0.18 + 0.65 * speed_norm(speed))
            ice_colors[..., :3] = 0.55 * ice_colors[..., :3] + 0.45 * np.array(
                [0.75, 0.9, 1.0]
            )
            ice_colors[..., 3] = np.where(has_ice, 0.9, 0.0)
            axis.plot_surface(x, y, ice_plot, facecolors=ice_colors,
                              linewidth=0.12, edgecolor="#d9f4ff",
                              antialiased=True, shade=True)

            arrow_stride = max(2, 8 // stride)
            selection = has_ice[::arrow_stride, ::arrow_stride]
            arrow_x = x[::arrow_stride, ::arrow_stride][selection]
            arrow_y = y[::arrow_stride, ::arrow_stride][selection]
            arrow_z = ice_z[::arrow_stride, ::arrow_stride][selection] + 0.12
            speed_sample = speed[::arrow_stride, ::arrow_stride][selection]
            axis.quiver(arrow_x, arrow_y, arrow_z,
                        np.ones_like(arrow_x), np.zeros_like(arrow_y),
                        np.zeros_like(arrow_z),
                        length=5.0, normalize=True, color="#d8f4ff",
                        alpha=np.clip(speed_sample / 60.0, 0.15, 0.75),
                        linewidth=0.65, arrow_length_ratio=0.3)

        # Highlight active glacial erosion without hiding the terrain.
        erosion = frame["glacial_erosion_m"]
        active = erosion > 0.15 * erosion_max
        if np.any(active):
            axis.scatter(x[active], y[active], terrain_z[active] + 0.04,
                         c=erosion[active], cmap="inferno", vmin=0,
                         vmax=erosion_max, s=5, alpha=0.65, depthshade=False)

        water_x, water_y = np.meshgrid([x.min(), x.max()], [y.min(), y.max()])
        axis.plot_surface(water_x, water_y, np.zeros_like(water_x),
                          color="#1b6d91", alpha=0.28, shade=False)

        axis.view_init(elev=31, azim=-62)
        axis.set_xlim(float(x.min()), float(x.max()))
        axis.set_ylim(float(y.min()), float(y.max()))
        axis.set_zlim(-6.0, 21.0)
        axis.set_box_aspect((1.6, 1.0, 0.48))
        axis.set_xlabel("distance toward coast (km)", color="#d9edf7", labelpad=10)
        axis.set_ylabel("across fjord (km)", color="#d9edf7", labelpad=10)
        axis.tick_params(colors="#9eb6c5", labelsize=8)
        axis.grid(False)
        for pane in (axis.xaxis.pane, axis.yaxis.pane, axis.zaxis.pane):
            pane.set_alpha(0.0)

        unloading = time_years > 100000.0
        phase = "POSTGLACIAL UNLOADING" if unloading else "GLACIAL INCISION"
        figure.suptitle("ICE · LANDSCAPE · SOLID-EARTH RESPONSE",
                        color="#f3f7fa", fontsize=20, fontweight="bold", y=0.97)
        status_axis.text(0.08, 0.93, phase, color="#64d5ff", fontsize=11,
                         fontweight="bold", transform=status_axis.transAxes)
        status_axis.text(0.08, 0.83, f"{time_years / 1000.0:5.1f} kyr",
                         color="#ffffff", fontsize=28, fontweight="bold",
                         transform=status_axis.transAxes)
        status_axis.text(0.08, 0.75, "model time", color="#9eb6c5", fontsize=10,
                         transform=status_axis.transAxes)

        metrics = (
            ("ICE THICKNESS", f"{np.max(frame['ice_thickness_m']):.0f} m"),
            ("BASAL SPEED", f"{np.max(frame['basal_ice_velocity_m_per_year']):.1f} m/yr"),
            ("EROSION / 100 yr", f"{np.max(erosion):.3f} m"),
            ("SEDIMENT FLUX", f"{np.max(frame['sediment_flux_m3_per_year']) / 1e6:.2f} Mm³/yr"),
            ("REBOUND SPEED", f"{np.max(frame['regional_ice_load_velocity_m_per_year']) * 1000:.1f} mm/yr"),
        )
        for row, (label, value) in enumerate(metrics):
            y_position = 0.62 - 0.105 * row
            status_axis.text(0.08, y_position, value, color="#f3f7fa", fontsize=16,
                             fontweight="bold", transform=status_axis.transAxes)
            status_axis.text(0.08, y_position - 0.035, label, color="#9eb6c5",
                             fontsize=8, transform=status_axis.transAxes)

        progress = time_years / 110000.0
        status_axis.plot([0.08, 0.92], [0.08, 0.08], color="#294456", linewidth=5,
                         solid_capstyle="round", transform=status_axis.transAxes)
        status_axis.plot([0.08, 0.08 + 0.84 * progress], [0.08, 0.08],
                         color="#64d5ff", linewidth=5, solid_capstyle="round",
                         transform=status_axis.transAxes)
        status_axis.text(0.08, 0.025, f"vertical exaggeration {VERTICAL_EXAGGERATION:.0f}×",
                         color="#7893a3", fontsize=8, transform=status_axis.transAxes)

    poster_frame = min(
        range(len(movie_frames)),
        key=lambda index: abs(movie_frames[index][1] - 50000.0),
    )
    draw(poster_frame)
    figure.savefig(poster, dpi=160, facecolor=figure.get_facecolor())

    if shutil.which("ffmpeg") is None:
        raise RuntimeError("ffmpeg is required to write the MP4 animation")
    writer = animation.FFMpegWriter(
        fps=frames_per_second,
        codec="libx264",
        bitrate=5000,
        extra_args=["-pix_fmt", "yuv420p"],
    )
    movie = animation.FuncAnimation(
        figure, draw, frames=len(movie_frames), interval=1000 / frames_per_second
    )
    movie.save(output, writer=writer, dpi=110)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stride", type=int, default=2,
                        help="render every Nth surface point")
    parser.add_argument("--frames-per-second", type=int, default=12)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "norway_fjord_cinematic.mp4")
    parser.add_argument("--poster", type=Path,
                        default=ROOT / "norway_fjord_cinematic.png")
    arguments = parser.parse_args()
    render(arguments.output, arguments.poster,
           max(1, arguments.stride), arguments.frames_per_second)


if __name__ == "__main__":
    main()
