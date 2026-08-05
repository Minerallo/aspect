#!/usr/bin/env python3
"""Analyze the local Fortran, C++, spherical, and climate verification runs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ASPECT_ROOT = HERE.parents[1]
WORKSPACE = HERE.parents[3]


def amplitude(values: np.ndarray, basis: np.ndarray) -> float:
    return float(np.dot(values, basis) / np.dot(basis, basis))


def relative_difference(first: float, second: float) -> float:
    return abs(first - second) / max(abs(first), abs(second), 1.0)


def read_last_surface(directory: Path) -> pd.DataFrame:
    files = sorted(directory.glob("surface-*.csv"))
    if not files:
        raise FileNotFoundError(f"No surface CSV files in {directory}")
    return pd.read_csv(files[-1])


def reshape_regular(data: pd.DataFrame, field: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.sort(data.x_m.unique())
    y = np.sort(data.y_m.unique())
    values = data.pivot(index="y_m", columns="x_m", values=field).loc[y, x].to_numpy()
    return x / 1000.0, y / 1000.0, values


def spherical_amplitudes(root: Path, dimension: int) -> tuple[float, float, float]:
    case = root / f"output-hill-diffusion-spherical-{dimension}d" / "solution"
    measured = []
    for output_index in (0, 2):
        data = np.loadtxt(case / f"solution_surface-{output_index:05d}.0000.gnuplot")
        xyz = data[:, :dimension]
        radius = np.linalg.norm(xyz, axis=1)
        outer = radius > 0.75
        xyz = xyz[outer]
        radius = radius[outer]
        displacement = radius - 1.0
        if dimension == 2:
            basis = np.cos(2.0 * np.arctan2(xyz[:, 1], xyz[:, 0]))
            eigenvalue = 4.0
        else:
            basis = 0.5 * (3.0 * (xyz[:, 2] / radius) ** 2 - 1.0)
            eigenvalue = 6.0
        measured.append(amplitude(displacement, basis))
    exact_final = 0.05 * math.exp(-0.1 * eigenvalue * 0.02)
    return measured[0], measured[1], exact_final


def rotation_metrics(directory: Path) -> dict[str, float]:
    surface = read_last_surface(directory / "fastscape_surface_evolution")
    longitude = np.deg2rad(surface.longitude_deg.to_numpy())
    latitude = np.deg2rad(surface.latitude_deg.to_numpy())
    expected = 1500.0 * np.sin(longitude) * np.sin(latitude)
    observed = surface.elevation_m.to_numpy()
    return {
        "cells": int(len(surface)),
        "correlation": float(np.corrcoef(expected, observed)[0, 1]),
        "amplitude_retained": float(
            np.sqrt(np.mean(observed**2)) / np.sqrt(np.mean(expected**2))
        ),
        "mean_elevation_m": float(np.mean(observed)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-tests", type=Path, default=HERE.parents[2] / "tests")
    parser.add_argument(
        "--spherical-tests",
        type=Path,
        default=WORKSPACE / "aspect_sphericaldiffusion" / "tests",
    )
    parser.add_argument(
        "--climate-loop",
        type=Path,
        default=ASPECT_ROOT / "cookbooks" / "fastscape_climberx_physical_loop"
        / "output-diagnostic-ice-sequence",
    )
    args = parser.parse_args()
    run_output = HERE / "output"
    results = HERE / "results"
    run_output.mkdir(exist_ok=True)
    results.mkdir(exist_ok=True)

    metrics: dict[str, object] = {}
    checks: dict[str, bool] = {}

    # Fortran analytical and process checks.
    diffusion_initial = pd.read_csv(run_output / "fortran_diffusion_initial.csv")
    diffusion_final = pd.read_csv(run_output / "fortran_diffusion_final.csv")
    length = 100e3
    diffusion_basis = np.sin(np.pi * diffusion_final.x_m / length) * np.sin(
        np.pi * diffusion_final.y_m / length
    )
    measured_diffusion = amplitude(diffusion_final.elevation_m.to_numpy(), diffusion_basis)
    exact_diffusion = 1000.0 * math.exp(-2.0 * math.pi**2 * 100.0 * 1e6 / length**2)
    diffusion_error = abs(measured_diffusion - exact_diffusion) / exact_diffusion
    metrics["fortran_diffusion"] = {
        "measured_amplitude_m": measured_diffusion,
        "exact_amplitude_m": exact_diffusion,
        "relative_error": diffusion_error,
    }
    checks["Fortran analytical diffusion"] = diffusion_error < 0.01

    river = pd.read_csv(run_output / "fortran_river_deposition_final.csv")
    sediment = river.elevation_m - river.basement_m
    metrics["fortran_river_deposition"] = {
        "maximum_erosion_rate_m_per_year": float(river.erosion_rate_m_per_year.max()),
        "maximum_sediment_thickness_m": float(sediment.max()),
        "maximum_drainage_area_m2": float(river.drainage_area_m2.max()),
    }
    checks["Fortran river erosion and deposition"] = bool(
        river.erosion_rate_m_per_year.max() > 0 and sediment.max() > 0
    )

    glacial_fortran = pd.read_csv(run_output / "fortran_glacial_final.csv")
    fortran_glacial_rate = float(glacial_fortran.erosion_rate_m_per_year.max())
    metrics["fortran_glacial"] = {
        "maximum_rate_m_per_year": fortran_glacial_rate,
        "saturated_law_expected_m_per_year": 0.01,
    }
    checks["Fortran glacial erosion law"] = abs(fortran_glacial_rate - 0.01) < 1e-10

    rotation_initial = pd.read_csv(run_output / "fortran_rotation_initial.csv")
    rotation_final = pd.read_csv(run_output / "fortran_rotation_final.csv")
    rotation_correlation = float(
        np.corrcoef(rotation_initial.elevation_m, rotation_final.elevation_m)[0, 1]
    )
    rotation_mass_error = abs(
        rotation_final.elevation_m.sum() / rotation_initial.elevation_m.sum() - 1.0
    )
    rotation_rms = float(
        np.sqrt(np.mean((rotation_final.elevation_m - rotation_initial.elevation_m) ** 2))
        / np.sqrt(np.mean(rotation_initial.elevation_m**2))
    )
    metrics["fortran_full_rotation"] = {
        "correlation": rotation_correlation,
        "relative_integral_error": rotation_mass_error,
        "relative_root_mean_square_error": rotation_rms,
        "peak_amplitude_retained": float(
            rotation_final.elevation_m.max() / rotation_initial.elevation_m.max()
        ),
    }
    checks["Fortran full-turn advection"] = (
        rotation_correlation > 0.98 and rotation_mass_error < 1e-3
    )

    # C++ process checks and serial/parallel equivalence.
    test_root = args.external_tests
    cpp_glacial = read_last_surface(
        test_root / "output-fastscape-cpp-box-glacial-erosion" / "fastscape_surface_evolution"
    )
    cpp_glacial_increment = float(cpp_glacial.glacial_erosion_m.max())
    metrics["cpp_glacial"] = {
        "maximum_erosion_per_500_year_substep_m": cpp_glacial_increment,
        "expected_m": 5.0,
    }
    checks["C++ glacial erosion law"] = abs(cpp_glacial_increment - 5.0) < 1e-10

    marine_dir = (
        test_root / "output-fastscape-cpp-global-marine-deposition"
        / "fastscape_surface_evolution"
    )
    cpp_marine = read_last_surface(marine_dir)
    marine_budget = pd.read_csv(marine_dir / "sediment_budget.csv")
    metrics["cpp_marine_deposition"] = {
        "maximum_sediment_thickness_m": float(cpp_marine.sediment_thickness_m.max()),
        "stored_solid_volume_m3": float(marine_budget.stored_sediment_solid_volume_m3.iloc[-1]),
    }
    checks["C++ marine deposition"] = bool(
        cpp_marine.sediment_thickness_m.max() > 0
        and marine_budget.stored_sediment_solid_volume_m3.iloc[-1] > 0
    )

    serial_budget = pd.read_csv(
        test_root / "output-fastscape-cpp-global" / "fastscape_surface_evolution"
        / "sediment_budget.csv"
    )
    parallel_budget = pd.read_csv(
        test_root / "output-fastscape-cpp-global-mpi" / "fastscape_surface_evolution"
        / "sediment_budget.csv"
    )
    compared_columns = [
        "eroded_volume_m3",
        "sediment_outflux_m3_per_year",
        "max_drainage_area_m2",
        "min_elevation_m",
        "max_elevation_m",
    ]
    parallel_error = max(
        relative_difference(float(serial_budget[column].iloc[-1]),
                            float(parallel_budget[column].iloc[-1]))
        for column in compared_columns
    )
    metrics["cpp_serial_parallel"] = {"maximum_relative_difference": parallel_error}
    checks["C++ serial/two-process equivalence"] = parallel_error < 1e-12

    zero_budget = pd.read_csv(
        run_output / "global-flat-surface-control"
        / "fastscape_surface_evolution" / "sediment_budget.csv"
    )
    checks["C++ flat-surface control"] = bool(
        np.all(zero_budget.eroded_volume_m3 == 0)
        and np.all(zero_budget.sediment_outflux_m3_per_year == 0)
    )

    stationary_budget = pd.read_csv(
        run_output / "global-rigid-rotation-control"
        / "fastscape_surface_evolution" / "sediment_budget.csv"
    )
    stationary_range = max(
        float(stationary_budget.min_elevation_m.max()
              - stationary_budget.min_elevation_m.min()),
        float(stationary_budget.max_elevation_m.max()
              - stationary_budget.max_elevation_m.min()),
    )
    metrics["cpp_stationary_surface"] = {
        "maximum_elevation_range_change_m": stationary_range
    }
    checks["C++ stationary-surface control"] = stationary_range < 1e-10

    rotation_cases = {
        "coarse": run_output / "global-rigid-rotation",
        "refined": run_output / "global-rigid-rotation-refined",
        "fine": run_output / "global-rigid-rotation-fine",
    }
    cpp_rotations = {name: rotation_metrics(path) for name, path in rotation_cases.items()}
    metrics["cpp_full_rotation"] = cpp_rotations
    checks["C++ full-turn mean preservation"] = all(
        abs(result["mean_elevation_m"]) < 0.1 for result in cpp_rotations.values()
    )
    checks["C++ full-turn convergence"] = all(
        cpp_rotations[name]["amplitude_retained"]
        < cpp_rotations[next_name]["amplitude_retained"]
        for name, next_name in (("coarse", "refined"), ("refined", "fine"))
    )

    # Independent spherical diffusion eigenfunction checks.
    spherical_results = {}
    for dimension in (2, 3):
        initial, final, exact = spherical_amplitudes(args.spherical_tests, dimension)
        error = abs(final - exact) / exact
        spherical_results[f"{dimension}d"] = {
            "initial_amplitude": initial,
            "final_amplitude": final,
            "exact_final_amplitude": exact,
            "relative_error": error,
        }
        checks[f"ASPECT spherical diffusion {dimension}D"] = error < 0.002
    metrics["spherical_diffusion"] = spherical_results

    # The completed lightweight CLIMBER-X loop is intentionally reused: a
    # climate year costs roughly 100 seconds on this machine.
    timing = json.loads((args.climate_loop / "timing-summary.json").read_text())
    final_climate_surface = Path(timing["final_surface"])
    metrics["climber_x_coupling"] = timing
    checks["CLIMBER-X diagnostic-ice coupling"] = bool(
        timing["coupling_windows"] >= 2
        and timing["climate_model_years"] >= 3
        and final_climate_surface.is_file()
    )

    checks = {name: bool(value) for name, value in checks.items()}
    metrics["checks"] = checks
    metrics["all_required_checks_pass"] = all(checks.values())
    (results / "verification-metrics.json").write_text(
        json.dumps(metrics, indent=2, default=str) + "\n", encoding="utf-8"
    )

    # Process-field figure.
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    process_panels = [
        (river, "elevation_m", "Fortran elevation after river erosion", "terrain"),
        (river.assign(sediment_thickness_m=sediment), "sediment_thickness_m",
         "Fortran deposited sediment", "cividis"),
        (glacial_fortran, "erosion_rate_m_per_year", "Fortran glacial erosion rate", "Blues"),
    ]
    for axis, (data, field, title, color_map) in zip(axes[0], process_panels):
        x, y, values = reshape_regular(data, field)
        image = axis.pcolormesh(x, y, values, shading="auto", cmap=color_map)
        axis.set(title=title, xlabel="x (km)", ylabel="y (km)")
        fig.colorbar(image, ax=axis, shrink=0.82)

    cpp_box = read_last_surface(
        test_root / "output-fastscape-cpp-box" / "fastscape_surface_evolution"
    )
    panels = [
        (cpp_box, "fluvial_erosion_m", "C++ river incision per step", "magma"),
        (cpp_glacial, "glacial_erosion_m", "C++ glacial incision per substep", "Blues"),
        (cpp_marine, "sediment_thickness_m", "C++ marine sediment thickness", "cividis"),
    ]
    for axis, (data, field, title, color_map) in zip(axes[1], panels):
        if title.startswith("C++ marine"):
            image = axis.scatter(data.longitude_deg, data.latitude_deg, c=data[field],
                                 s=14, cmap=color_map)
            axis.set(xlabel="longitude (degrees)", ylabel="latitude (degrees)")
        else:
            image = axis.scatter(data.surface_x_m / 1000.0, data.surface_y_m / 1000.0,
                                 c=data[field], s=22, cmap=color_map)
            axis.set(xlabel="x (km)", ylabel="y (km)")
        axis.set_title(title)
        fig.colorbar(image, ax=axis, shrink=0.82)
    fig.suptitle("FastScape process verification: Fortran and C++")
    fig.savefig(results / "cross-code-processes.png", dpi=180)
    plt.close(fig)

    # Numerical behavior figure.
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    x, y, initial_grid = reshape_regular(rotation_initial, "elevation_m")
    _, _, final_grid = reshape_regular(rotation_final, "elevation_m")
    center = len(y) // 2
    axes[0, 0].plot(x, initial_grid[center], label="initial")
    axes[0, 0].plot(x, final_grid[center], label="after one turn")
    axes[0, 0].set(title="Fortran full rotation", xlabel="x (km)", ylabel="elevation (m)")
    axes[0, 0].legend()

    names = list(cpp_rotations)
    cells = [cpp_rotations[name]["cells"] for name in names]
    retained = [100.0 * cpp_rotations[name]["amplitude_retained"] for name in names]
    axes[0, 1].plot(cells, retained, "o-")
    axes[0, 1].set_xscale("log", base=2)
    axes[0, 1].set(title="C++ spherical advection convergence",
                   xlabel="surface cells", ylabel="amplitude retained after one turn (%)")
    for count, value in zip(cells, retained):
        axes[0, 1].annotate(f"{value:.1f}%", (count, value), xytext=(0, 7),
                            textcoords="offset points", ha="center")

    dimensions = ["2D", "3D"]
    exact_values = [spherical_results["2d"]["exact_final_amplitude"],
                    spherical_results["3d"]["exact_final_amplitude"]]
    observed_values = [spherical_results["2d"]["final_amplitude"],
                       spherical_results["3d"]["final_amplitude"]]
    positions = np.arange(2)
    axes[1, 0].bar(positions - 0.18, exact_values, 0.36, label="analytical")
    axes[1, 0].bar(positions + 0.18, observed_values, 0.36, label="ASPECT")
    axes[1, 0].set_xticks(positions, dimensions)
    axes[1, 0].set(title="Spherical diffusion eigenfunction decay", ylabel="final amplitude")
    axes[1, 0].legend()

    labels = list(checks)
    values = [1 if checks[label] else 0 for label in labels]
    colors = ["#2a9d8f" if value else "#e76f51" for value in values]
    axes[1, 1].barh(np.arange(len(labels)), values, color=colors)
    axes[1, 1].set_yticks(np.arange(len(labels)), labels)
    axes[1, 1].set_xlim(0, 1.05)
    axes[1, 1].set_xticks([0, 1], ["fail", "pass"])
    axes[1, 1].set_title("Required verification checks")
    fig.suptitle("Advection, diffusion, conservation, and coupled-model checks")
    fig.savefig(results / "numerical-verification.png", dpi=180)
    plt.close(fig)

    print(json.dumps({"all_pass": all(checks.values()), "checks": checks}, indent=2))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
