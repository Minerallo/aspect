#!/usr/bin/env python3
"""Extract and plot the paired local uplift experiment."""

import argparse
from pathlib import Path


FIELDS = ("elevation", "sediment_flux", "glacial_erosion", "ice_thickness")


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
    for cell_index in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(cell_index)
        point_ids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
        centers[cell_index] = points[point_ids, :2].mean(axis=0)

    arrays = {field: np.empty(grid.GetNumberOfCells()) for field in FIELDS}
    for field in FIELDS:
        array = grid.GetPointData().GetArray(field)
        if array is None:
            raise RuntimeError(f"{field!r} is missing from {path}")
        point_values = vtk_to_numpy(array)
        for cell_index in range(grid.GetNumberOfCells()):
            cell = grid.GetCell(cell_index)
            point_ids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
            arrays[field][cell_index] = point_values[point_ids].mean()
    return centers, arrays


def extract(without_glacial: Path, with_glacial: Path, output: Path) -> None:
    import numpy as np

    centers_without, fields_without = read_cell_data(without_glacial)
    centers_with, fields_with = read_cell_data(with_glacial)
    np.testing.assert_allclose(centers_without, centers_with)
    output.parent.mkdir(parents=True, exist_ok=True)
    values = {"centers": centers_without}
    for field in FIELDS:
        values[f"without_{field}"] = fields_without[field]
        values[f"with_{field}"] = fields_with[field]
    np.savez_compressed(output, **values)


def plot(input_path: Path, output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    import numpy as np
    from matplotlib.colors import TwoSlopeNorm

    data = np.load(input_path)
    centers = data["centers"] / 1000.0
    triangulation = mtri.Triangulation(centers[:, 0], centers[:, 1])
    elevation_without = data["without_elevation"]
    elevation_with = data["with_elevation"]
    elevation_difference = elevation_with - elevation_without
    flux_without = np.log10(1.0 + data["without_sediment_flux"])
    flux_with = np.log10(1.0 + data["with_sediment_flux"])
    glacial_erosion = data["with_glacial_erosion"]

    figure, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    elevation_limits = (min(elevation_without.min(), elevation_with.min()),
                        max(elevation_without.max(), elevation_with.max()))
    flux_limits = (min(flux_without.min(), flux_with.min()),
                   max(flux_without.max(), flux_with.max()))
    panels = (
        (elevation_without, "Without glacial erosion: elevation (m)",
         "terrain", elevation_limits, None),
        (elevation_with, "With glacial erosion: elevation (m)",
         "terrain", elevation_limits, None),
        (elevation_difference, "Elevation difference: with minus without (m)",
         "coolwarm", None,
         TwoSlopeNorm(vcenter=0.0,
                      vmin=elevation_difference.min(),
                      vmax=max(1.e-12, elevation_difference.max()))),
        (flux_without, "Without glacial erosion: log₁₀(1 + sediment flux)",
         "viridis", flux_limits, None),
        (flux_with, "With glacial erosion: log₁₀(1 + sediment flux)",
         "viridis", flux_limits, None),
        (glacial_erosion, "Glacial erosion in final 500-year landscape step (m)",
         "magma", None, None),
    )
    for axis, (values, title, color_map, limits, normalization) in zip(axes.flat, panels):
        options = {"levels": 24, "cmap": color_map, "extend": "both"}
        if limits is not None:
            options["levels"] = np.linspace(limits[0], limits[1], 24)
        if normalization is not None:
            options["norm"] = normalization
        surface = axis.tricontourf(triangulation, values, **options)
        axis.tricontour(triangulation, data["with_ice_thickness"],
                        levels=[500.0], colors="cyan", linewidths=1.0)
        axis.set(title=title, xlabel="x (km)", ylabel="y (km)", aspect="equal")
        figure.colorbar(surface, ax=axis, shrink=0.82)

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    extract_parser = subparsers.add_parser("extract")
    extract_parser.add_argument("without_glacial", type=Path)
    extract_parser.add_argument("with_glacial", type=Path)
    extract_parser.add_argument("output", type=Path)
    plot_parser = subparsers.add_parser("plot")
    plot_parser.add_argument("input", type=Path)
    plot_parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    if arguments.command == "extract":
        extract(arguments.without_glacial, arguments.with_glacial, arguments.output)
    else:
        plot(arguments.input, arguments.output)


if __name__ == "__main__":
    main()
