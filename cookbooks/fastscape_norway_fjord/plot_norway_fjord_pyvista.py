#!/usr/bin/env python3
"""Render the Norway-like glaciation and unloading experiment with PyVista."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pyvista as pv


ROOT = Path(__file__).resolve().parent


def read(output: str, step: int) -> dict[str, np.ndarray]:
    path = ROOT / output / "fastscape_surface_evolution" / f"surface-{step:05d}.csv"
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {
        field: np.asarray([float(row[field]) for row in rows])
        for field in rows[0]
        if field != "bedrock_lithology"
    }


def surface(data: dict[str, np.ndarray], height: np.ndarray, scalar_name: str,
            scalars: np.ndarray) -> pv.PolyData:
    points = np.column_stack(
        (data["surface_x_m"] / 1000.0, data["surface_y_m"] / 1000.0, height / 1000.0)
    )
    mesh = pv.PolyData(points).delaunay_2d()
    mesh[scalar_name] = scalars
    return mesh


loaded = read("output-fastscape-norway-glacial", 20)
unloaded = read("output-fastscape-norway-glacial", 22)
control = read("output-fastscape-norway-control", 20)
incision = control["elevation_m"] - loaded["elevation_m"]
rebound = unloaded["elevation_m"] - loaded["elevation_m"]

terrain = surface(loaded, loaded["elevation_m"], "elevation (m)", loaded["elevation_m"])
ice_height = loaded["elevation_m"] + loaded["ice_thickness_m"]
ice_mesh = surface(loaded, ice_height, "ice thickness (m)", loaded["ice_thickness_m"])
ice_mesh = (
    ice_mesh.threshold(100.0, scalars="ice thickness (m)")
    .extract_surface(algorithm="dataset_surface")
    .triangulate()
)
erosion_mesh = surface(loaded, loaded["elevation_m"], "last-step erosion (m)", loaded["glacial_erosion_m"])
incision_mesh = surface(loaded, loaded["elevation_m"], "extra incision (m)", incision)
rebound_mesh = surface(unloaded, unloaded["elevation_m"], "rebound (m)", rebound)

# Preserve the actual PyVista meshes as a reusable interactive data product.
blocks = pv.MultiBlock(
    {
        "glaciated_terrain": terrain,
        "ice_surface": ice_mesh,
        "last_step_glacial_erosion": erosion_mesh,
        "incision_relative_to_control": incision_mesh,
        "unloaded_rebound": rebound_mesh,
    }
)
blocks.save(ROOT / "norway_fjord_results.vtm")


def native_pyvista_render() -> None:
    """Use VTK rendering where a graphical or OSMesa backend is available."""
    pv.global_theme.font.size = 11
    pv.global_theme.font.title_size = 13
    plotter = pv.Plotter(shape=(2, 2), off_screen=True, window_size=(1800, 1100))
    camera = [(225.0, -175.0, 105.0), (80.0, 50.0, 0.0), (0.0, 0.0, 1.0)]
    panels = (
        (terrain, "elevation (m)", "terrain", "Glaciated fjord after 100,000 years"),
        (erosion_mesh, "last-step erosion (m)", "inferno", "Glacial erosion in the last 100-year step"),
        (incision_mesh, "extra incision (m)", "magma", "Difference from no-glacial-erosion control"),
        (rebound_mesh, "rebound (m)", "viridis", "Landscape after 10,000 years of unloading"),
    )
    for index, (mesh, scalars, colors, title) in enumerate(panels):
        plotter.subplot(index // 2, index % 2)
        plotter.add_mesh(mesh, scalars=scalars, cmap=colors)
        if index == 0:
            plotter.add_mesh(ice_mesh, color="#bfe8ff", opacity=0.58,
                             show_scalar_bar=False)
        plotter.add_text(title, position="upper_left")
        plotter.camera_position = camera
        plotter.show_axes()
    plotter.link_views()
    plotter.screenshot(ROOT / "norway_fjord_pyvista.png")
    plotter.close()


def headless_render() -> None:
    """Render PyVista-prepared triangles when native VTK lacks a display."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(16, 9.5), constrained_layout=True)
    panels = (
        (terrain, "elevation (m)", "terrain", "Glaciated fjord after 100,000 years"),
        (erosion_mesh, "last-step erosion (m)", "inferno", "Glacial erosion in the last 100-year step"),
        (incision_mesh, "extra incision (m)", "magma", "Difference from no-glacial-erosion control"),
        (rebound_mesh, "rebound (m)", "viridis", "Landscape after 10,000 years of unloading"),
    )
    for index, (mesh, scalar, colors, title) in enumerate(panels, start=1):
        axis = figure.add_subplot(2, 2, index, projection="3d")
        triangles = mesh.faces.reshape((-1, 4))[:, 1:]
        artist = axis.plot_trisurf(
            mesh.points[:, 0], mesh.points[:, 1], mesh.points[:, 2],
            triangles=triangles, cmap=colors, linewidth=0,
            antialiased=True, shade=False,
        )
        face_values = mesh[scalar][triangles].mean(axis=1)
        artist.set_array(face_values)
        artist.set_clim(float(np.min(face_values)), float(np.max(face_values)))
        if index == 1:
            ice_triangles = ice_mesh.faces.reshape((-1, 4))[:, 1:]
            axis.plot_trisurf(
                ice_mesh.points[:, 0], ice_mesh.points[:, 1], ice_mesh.points[:, 2],
                triangles=ice_triangles, color="#bfe8ff", alpha=0.5,
                linewidth=0, shade=True,
            )
        # A translucent sea-level plane makes the overdeepened outlet legible.
        xx, yy = np.meshgrid([0.0, 160.0], [0.0, 100.0])
        axis.plot_surface(xx, yy, np.zeros_like(xx), color="#65bde6", alpha=0.12)
        figure.colorbar(artist, ax=axis, shrink=0.58, pad=0.09, label=scalar)
        axis.set_title(title)
        axis.set_xlabel("x (km)")
        axis.set_ylabel("y (km)")
        axis.view_init(elev=28, azim=-58)
        axis.set_box_aspect((1.6, 1.0, 0.35))
    figure.suptitle(
        "Synthetic western Norway process test — meshes constructed and exported with PyVista"
    )
    figure.savefig(ROOT / "norway_fjord_pyvista.png", dpi=180)
    plt.close(figure)


if "--native-pyvista" in sys.argv:
    native_pyvista_render()
else:
    headless_render()
