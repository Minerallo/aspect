#!/usr/bin/env python3
"""Plot the FastScape C++ lithology and provenance demonstration."""

import argparse
from pathlib import Path


def read_surface(path: Path):
    import numpy as np
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    points = vtk_to_numpy(grid.GetPoints().GetData())
    field_names = (
        "elevation",
        "bedrock_lithology",
        "sediment_thickness",
        "sediment_fraction_limestone",
    )
    point_fields = {
        name: vtk_to_numpy(grid.GetPointData().GetArray(name))
        for name in field_names
    }
    centers = np.empty((grid.GetNumberOfCells(), 2))
    fields = {name: np.empty(grid.GetNumberOfCells()) for name in field_names}
    for cell_index in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(cell_index)
        point_ids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
        centers[cell_index] = points[point_ids, :2].mean(axis=0)
        for name in field_names:
            fields[name][cell_index] = point_fields[name][point_ids].mean()
    return centers, fields


def to_grid(centers, values):
    import numpy as np

    x = np.unique(centers[:, 0])
    y = np.unique(centers[:, 1])
    grid = np.empty((len(y), len(x)))
    grid[np.searchsorted(y, centers[:, 1]), np.searchsorted(x, centers[:, 0])] = values
    return x / 1000.0, y / 1000.0, grid


def plot(input_path: Path, output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import ListedColormap

    centers, fields = read_surface(input_path)
    panels = []
    for name in fields:
        x, y, values = to_grid(centers, fields[name])
        panels.append((name, values))
    extent = (x.min(), x.max(), y.min(), y.max())

    figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    definitions = (
        (panels[1][1], "Bedrock lithology (0 granite, 1 limestone)",
         ListedColormap(["#d98b5f", "#c9d8a4"]), 0.0, 1.0),
        (panels[0][1], "Elevation (m)", "terrain", None, None),
        (np.log10(np.maximum(panels[2][1], 1e-6)),
         "Log₁₀ sediment thickness (m)", "viridis", -6.0, None),
        (panels[3][1], "Limestone fraction in deposited sediment",
         "cividis", 0.0, 1.0),
    )
    for axis, (values, title, color_map, minimum, maximum) in zip(
            axes.flat, definitions):
        image = axis.imshow(values, origin="lower", extent=extent,
                            cmap=color_map, vmin=minimum, vmax=maximum,
                            interpolation="nearest", aspect="equal")
        axis.set(title=title, xlabel="x (km)", ylabel="y (km)")
        figure.colorbar(image, ax=axis, shrink=0.82)

    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    plot(arguments.input, arguments.output)


if __name__ == "__main__":
    main()
