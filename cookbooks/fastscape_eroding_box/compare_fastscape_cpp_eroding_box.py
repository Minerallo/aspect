#!/usr/bin/env python3
"""Extract and plot the paired FastScape C++ eroding-box runs."""

import argparse
from pathlib import Path


FIELDS = (
    "elevation",
    "drainage_area",
    "fluvial_erosion",
    "glacial_erosion",
    "sediment_flux",
    "ice_thickness",
)


def read_cell_data(path: Path):
    import numpy as np
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    points = vtk_to_numpy(grid.GetPoints().GetData())
    centers = np.empty((grid.GetNumberOfCells(), 2))
    arrays = {field: np.empty(grid.GetNumberOfCells()) for field in FIELDS}
    point_arrays = {
        field: vtk_to_numpy(grid.GetPointData().GetArray(field)) for field in FIELDS
    }
    for cell_index in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(cell_index)
        point_ids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
        centers[cell_index] = points[point_ids, :2].mean(axis=0)
        for field in FIELDS:
            arrays[field][cell_index] = point_arrays[field][point_ids].mean()
    return centers, arrays


def extract(normal_path: Path, glacial_path: Path, output: Path) -> None:
    import numpy as np

    normal_centers, normal_fields = read_cell_data(normal_path)
    glacial_centers, glacial_fields = read_cell_data(glacial_path)
    np.testing.assert_allclose(normal_centers, glacial_centers)
    values = {"centers": normal_centers}
    for field in FIELDS:
        values[f"normal_{field}"] = normal_fields[field]
        values[f"glacial_{field}"] = glacial_fields[field]
    np.savez_compressed(output, **values)


def to_grid(centers, values):
    import numpy as np

    x = np.unique(centers[:, 0])
    y = np.unique(centers[:, 1])
    grid = np.empty((len(y), len(x)))
    grid[np.searchsorted(y, centers[:, 1]), np.searchsorted(x, centers[:, 0])] = values
    return x / 1000.0, y / 1000.0, grid


def plot(input_path: Path, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import TwoSlopeNorm

    data = np.load(input_path)
    centers = data["centers"]
    x, y, normal_elevation = to_grid(centers, data["normal_elevation"])
    _, _, glacial_elevation = to_grid(centers, data["glacial_elevation"])
    _, _, normal_drainage = to_grid(
        centers, np.log10(np.maximum(data["normal_drainage_area"], 1.0))
    )
    _, _, glacial_drainage = to_grid(
        centers, np.log10(np.maximum(data["glacial_drainage_area"], 1.0))
    )
    _, _, glacial_erosion = to_grid(centers, data["glacial_glacial_erosion"])
    _, _, ice = to_grid(centers, data["glacial_ice_thickness"])
    difference = glacial_elevation - normal_elevation
    extent = (x.min(), x.max(), y.min(), y.max())

    figure, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    elevation_min = min(normal_elevation.min(), glacial_elevation.min())
    elevation_max = max(normal_elevation.max(), glacial_elevation.max())
    drainage_min = min(normal_drainage.min(), glacial_drainage.min())
    drainage_max = max(normal_drainage.max(), glacial_drainage.max())
    panels = (
        (normal_elevation, "River erosion: elevation (m)", "terrain",
         elevation_min, elevation_max, None),
        (glacial_elevation, "River and glacial erosion: elevation (m)", "terrain",
         elevation_min, elevation_max, None),
        (difference, "Elevation difference: glacial minus river-only (m)", "coolwarm",
         -150.0, 150.0, TwoSlopeNorm(vcenter=0.0, vmin=-150.0, vmax=150.0)),
        (normal_drainage, "River erosion: log₁₀ drainage area (m²)", "Blues",
         drainage_min, drainage_max, None),
        (glacial_drainage, "River and glacial erosion: log₁₀ drainage area (m²)", "Blues",
         drainage_min, drainage_max, None),
        (glacial_erosion, "Glacial erosion in final 10,000-year step (m)", "magma",
         glacial_erosion.min(), glacial_erosion.max(), None),
    )

    for index, (axis, panel) in enumerate(zip(axes.flat, panels)):
        values, title, color_map, minimum, maximum, normalization = panel
        image = axis.imshow(
            values,
            origin="lower",
            extent=extent,
            cmap=color_map,
            vmin=None if normalization else minimum,
            vmax=None if normalization else maximum,
            norm=normalization,
            interpolation="nearest",
            aspect="equal",
        )
        if index < 2:
            drainage = normal_drainage if index == 0 else glacial_drainage
            axis.contour(x, y, drainage, levels=[8.0, 8.5, 9.0],
                         colors=["#5bc0eb", "#168aad", "#023e8a"],
                         linewidths=[0.35, 0.7, 1.1])
        if index in (1, 2, 4, 5):
            axis.contour(x, y, ice, levels=[500.0], colors="cyan", linewidths=0.8)
        axis.set(title=title, xlabel="x (km)", ylabel="y (km)")
        figure.colorbar(image, ax=axis, shrink=0.82)

    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    extract_parser = subparsers.add_parser("extract")
    extract_parser.add_argument("normal", type=Path)
    extract_parser.add_argument("glacial", type=Path)
    extract_parser.add_argument("output", type=Path)
    plot_parser = subparsers.add_parser("plot")
    plot_parser.add_argument("input", type=Path)
    plot_parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    if arguments.command == "extract":
        extract(arguments.normal, arguments.glacial, arguments.output)
    else:
        plot(arguments.input, arguments.output)


if __name__ == "__main__":
    main()
