#!/usr/bin/env python3
"""Native surface exchange for CLIMBER-X, FastScape, and ASPECT.

The format is deliberately small and dependency-free: a fixed binary header,
longitude/latitude point coordinates, and named double-precision fields.  It
is intended for coupling state, not archival model output; NetCDF remains the
right format for diagnostics and long-term restarts.
"""

from __future__ import annotations

import argparse
import csv
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np


MAGIC = b"CXCHG001"
VERSION = 1
HEADER_LAYOUT = "8sIIIId"
FIELD_TEXT_LENGTH = 32
SECONDS_PER_YEAR = 365.0 * 24.0 * 3600.0


@dataclass
class SurfaceExchange:
    model_time_years: float
    longitude: np.ndarray
    latitude: np.ndarray
    fields: dict[str, np.ndarray]
    units: dict[str, str]
    byte_order: str = "<"


def _fixed_text(value: str) -> bytes:
    encoded = value.encode("ascii")
    if len(encoded) > FIELD_TEXT_LENGTH:
        raise ValueError(f"exchange metadata is longer than {FIELD_TEXT_LENGTH}: {value}")
    return encoded.ljust(FIELD_TEXT_LENGTH, b" ")


def write_exchange(path: Path, exchange: SurfaceExchange) -> None:
    longitude = np.asarray(exchange.longitude, dtype="<f8").reshape(-1)
    latitude = np.asarray(exchange.latitude, dtype="<f8").reshape(-1)
    if longitude.shape != latitude.shape:
        raise ValueError("longitude and latitude must contain the same number of points")
    number_of_points = longitude.size

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        header = struct.Struct(exchange.byte_order + HEADER_LAYOUT)
        stream.write(
            header.pack(
                MAGIC,
                VERSION,
                number_of_points,
                len(exchange.fields),
                0,
                exchange.model_time_years,
            )
        )
        numeric_dtype = np.dtype(exchange.byte_order + "f8")
        stream.write(longitude.astype(numeric_dtype, copy=False).tobytes())
        stream.write(latitude.astype(numeric_dtype, copy=False).tobytes())
        for name, field in exchange.fields.items():
            values = np.asarray(field, dtype="<f8").reshape(-1)
            if values.size != number_of_points:
                raise ValueError(f"{name} has {values.size} values; expected {number_of_points}")
            stream.write(_fixed_text(name))
            stream.write(_fixed_text(exchange.units.get(name, "")))
            stream.write(values.astype(numeric_dtype, copy=False).tobytes())
    temporary.replace(path)


def read_exchange(path: Path) -> SurfaceExchange:
    with path.open("rb") as stream:
        header_size = struct.calcsize(HEADER_LAYOUT)
        raw_header = stream.read(header_size)
        if len(raw_header) != header_size:
            raise ValueError(f"{path}: truncated surface-exchange header")
        little_values = struct.unpack("<" + HEADER_LAYOUT, raw_header)
        big_values = struct.unpack(">" + HEADER_LAYOUT, raw_header)
        if little_values[0] == MAGIC and little_values[1] == VERSION:
            byte_order = "<"
            magic, version, number_of_points, number_of_fields, _, model_time = little_values
        elif big_values[0] == MAGIC and big_values[1] == VERSION:
            byte_order = ">"
            magic, version, number_of_points, number_of_fields, _, model_time = big_values
        else:
            raise ValueError(f"{path}: unsupported surface-exchange header")

        longitude = np.fromfile(stream, dtype=byte_order + "f8", count=number_of_points)
        latitude = np.fromfile(stream, dtype=byte_order + "f8", count=number_of_points)
        fields: dict[str, np.ndarray] = {}
        units: dict[str, str] = {}
        for _ in range(number_of_fields):
            name = stream.read(FIELD_TEXT_LENGTH).decode("ascii").strip()
            unit = stream.read(FIELD_TEXT_LENGTH).decode("ascii").strip()
            values = np.fromfile(stream, dtype=byte_order + "f8", count=number_of_points)
            if values.size != number_of_points:
                raise ValueError(f"{path}: truncated field {name}")
            fields[name] = values
            units[name] = unit
        if stream.read(1):
            raise ValueError(f"{path}: unexpected bytes after the final field")

    return SurfaceExchange(model_time, longitude, latitude, fields, units, byte_order)


def climate_controls(exchange: SurfaceExchange) -> tuple[np.ndarray, np.ndarray]:
    precipitation = np.maximum(exchange.fields["precipitation_rate"], 0.0)
    temperature = exchange.fields["surface_temperature"]
    weights = np.maximum(np.cos(np.deg2rad(exchange.latitude)), 0.0)
    wet = precipitation > 0.0
    reference = np.sum(precipitation[wet] * weights[wet]) / np.sum(weights[wet])
    runoff = np.clip(precipitation / reference, 0.0, 5.0)
    temperature_factor = np.clip(2.0 ** ((temperature - 288.15) / 10.0), 0.1, 3.0)
    erosion_strength = np.clip(np.maximum(runoff, 0.05) ** 0.4 * temperature_factor, 0.02, 5.0)
    return erosion_strength, runoff


def _regular_grid(exchange: SurfaceExchange) -> tuple[np.ndarray, np.ndarray]:
    first_latitude = exchange.latitude[0]
    changes = np.flatnonzero(np.abs(exchange.latitude - first_latitude) > 1.0e-10)
    nx = int(changes[0]) if changes.size else exchange.latitude.size
    if nx < 2 or exchange.longitude.size % nx:
        raise ValueError("exchange points do not form a regular longitude-latitude grid")
    ny = exchange.longitude.size // nx
    return exchange.longitude[:nx], exchange.latitude[::nx][:ny]


def write_aspect_structured(path: Path, exchange: SurfaceExchange, values: np.ndarray) -> None:
    longitude, latitude = _regular_grid(exchange)
    source_longitude = np.mod(exchange.longitude, 360.0)
    source_colatitude = 90.0 - exchange.latitude
    order = np.lexsort((source_longitude, source_colatitude))
    grid_values = np.asarray(values)[order].reshape(latitude.size, longitude.size)

    # ASPECT's spherical structured-data coordinates are longitude and
    # colatitude in radians. Add a periodic seam and pole values so its
    # interpolator sees a complete closed spherical domain.
    seam = 0.5 * (grid_values[:, 0] + grid_values[:, -1])
    values_with_seam = np.column_stack((seam, grid_values, seam))
    north_pole = np.full(values_with_seam.shape[1], np.mean(values_with_seam[0]))
    south_pole = np.full(values_with_seam.shape[1], np.mean(values_with_seam[-1]))
    closed_values = np.vstack((north_pole, values_with_seam, south_pole))
    closed_longitude = np.concatenate(([0.0], np.unique(source_longitude), [360.0]))
    closed_colatitude = np.concatenate(([0.0], np.unique(source_colatitude), [180.0]))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(f"# POINTS: {closed_longitude.size} {closed_colatitude.size}\n")
        stream.write("longitude colatitude value\n")
        for j, colatitude in enumerate(closed_colatitude):
            for i, lon in enumerate(closed_longitude):
                stream.write(
                    f"{np.deg2rad(lon):.16g} {np.deg2rad(colatitude):.16g} "
                    f"{closed_values[j, i]:.16g}\n"
                )


def spin_axis_from_history(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{path}: polar-wander history is empty")
    longitude = np.deg2rad(float(rows[-1]["pole_longitude_degrees"]))
    latitude = np.deg2rad(float(rows[-1]["pole_latitude_degrees"]))
    return np.array(
        [
            np.cos(latitude) * np.cos(longitude),
            np.cos(latitude) * np.sin(longitude),
            np.sin(latitude),
        ]
    )


def climate_vectors_in_body_frame(
    longitude_degrees: np.ndarray,
    latitude_degrees: np.ndarray,
    spin_axis: np.ndarray,
) -> np.ndarray:
    climate_north = np.asarray(spin_axis, dtype=float)
    climate_north /= np.linalg.norm(climate_north)
    climate_zero_longitude = np.array([1.0, 0.0, 0.0])
    climate_zero_longitude -= (
        climate_zero_longitude @ climate_north
    ) * climate_north
    if np.linalg.norm(climate_zero_longitude) < 1.0e-12:
        climate_zero_longitude = np.array([0.0, 1.0, 0.0])
        climate_zero_longitude -= (
            climate_zero_longitude @ climate_north
        ) * climate_north
    climate_zero_longitude /= np.linalg.norm(climate_zero_longitude)
    climate_east = np.cross(climate_north, climate_zero_longitude)

    longitude = np.deg2rad(longitude_degrees)
    latitude = np.deg2rad(latitude_degrees)
    return (
        (np.cos(latitude) * np.cos(longitude))[:, None] * climate_zero_longitude
        + (np.cos(latitude) * np.sin(longitude))[:, None] * climate_east
        + np.sin(latitude)[:, None] * climate_north
    )


def interpolate_spherical_surface(
    source_vectors: np.ndarray,
    source_values: np.ndarray,
    target_vectors: np.ndarray,
    number_of_neighbors: int,
) -> np.ndarray:
    if number_of_neighbors < 1:
        raise ValueError("the number of interpolation neighbors must be positive")
    number_of_neighbors = min(number_of_neighbors, source_vectors.shape[0])
    result = np.empty(target_vectors.shape[0])
    for start in range(0, target_vectors.shape[0], 256):
        stop = min(start + 256, target_vectors.shape[0])
        cosine = np.clip(target_vectors[start:stop] @ source_vectors.T, -1.0, 1.0)
        neighbors = np.argpartition(
            -cosine, number_of_neighbors-1, axis=1
        )[:, :number_of_neighbors]
        neighbor_cosine = np.take_along_axis(cosine, neighbors, axis=1)
        chord_distance = np.sqrt(np.maximum(2.0-2.0*neighbor_cosine, 0.0))
        exact = chord_distance < 1.0e-12
        weights = 1.0 / np.maximum(chord_distance, 1.0e-12)
        if np.any(exact):
            weights[exact.any(axis=1)] = exact[exact.any(axis=1)]
        weights /= weights.sum(axis=1, keepdims=True)
        result[start:stop] = np.sum(source_values[neighbors] * weights, axis=1)
    return result


def surface_to_climate_exchange(
    surface_csv: Path,
    climate_file: Path,
    output: Path,
    spin_axis: np.ndarray | None = None,
    number_of_neighbors: int = 4,
    reference_surface_csv: Path | None = None,
) -> None:
    climate = read_exchange(climate_file)
    with surface_csv.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    surface_longitude = np.deg2rad(np.array([float(row["longitude_deg"]) for row in rows]))
    surface_latitude = np.deg2rad(np.array([float(row["latitude_deg"]) for row in rows]))
    surface_change = np.array([float(row["elevation_change_m"]) for row in rows])
    if reference_surface_csv is not None:
        with reference_surface_csv.open(newline="", encoding="utf-8") as stream:
            reference_rows = list(csv.DictReader(stream))
        if len(reference_rows) != len(rows):
            raise ValueError("current and reference surfaces contain different cell counts")
        reference_longitude = np.array(
            [float(row["longitude_deg"]) for row in reference_rows]
        )
        reference_latitude = np.array(
            [float(row["latitude_deg"]) for row in reference_rows]
        )
        if not (
            np.allclose(np.rad2deg(surface_longitude), reference_longitude)
            and np.allclose(np.rad2deg(surface_latitude), reference_latitude)
        ):
            raise ValueError("current and reference surface cells do not match")
        surface_change -= np.array(
            [float(row["elevation_change_m"]) for row in reference_rows]
        )

    surface_vectors = np.column_stack(
        (
            np.cos(surface_latitude) * np.cos(surface_longitude),
            np.cos(surface_latitude) * np.sin(surface_longitude),
            np.sin(surface_latitude),
        )
    )
    if spin_axis is None:
        spin_axis = np.array([0.0, 0.0, 1.0])
    climate_vectors = climate_vectors_in_body_frame(
        climate.longitude,
        climate.latitude,
        spin_axis,
    )
    remapped_surface_change = interpolate_spherical_surface(
        surface_vectors,
        surface_change,
        climate_vectors,
        number_of_neighbors,
    )

    write_exchange(
        output,
        SurfaceExchange(
            climate.model_time_years,
            climate.longitude,
            climate.latitude,
            {"elevation_change": remapped_surface_change},
            {"elevation_change": "m"},
            climate.byte_order,
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    inspect_parser = subparsers.add_parser("inspect")
    inspect_parser.add_argument("exchange", type=Path)
    climate_parser = subparsers.add_parser("climate-to-aspect")
    climate_parser.add_argument("exchange", type=Path)
    climate_parser.add_argument(
        "--surface-topography",
        type=Path,
        help="write CLIMBER-X surface elevation as ASPECT spherical input",
    )
    climate_parser.add_argument("--erosion-strength", type=Path, required=True)
    climate_parser.add_argument("--surface-runoff", type=Path, required=True)
    climate_parser.add_argument("--ice-thickness", type=Path)
    climate_parser.add_argument("--basal-ice-velocity", type=Path)
    climate_parser.add_argument(
        "--prescribed-basal-ice-velocity",
        type=float,
        help=(
            "replace unavailable dynamic basal velocity with this explicit "
            "constant in grounded-ice cells, in meters per year"
        ),
    )
    return_parser = subparsers.add_parser("surface-to-climate")
    return_parser.add_argument("--surface", type=Path, required=True)
    return_parser.add_argument("--climate", type=Path, required=True)
    return_parser.add_argument("--output", type=Path, required=True)
    return_parser.add_argument(
        "--reference-surface",
        type=Path,
        help=(
            "subtract the cumulative elevation change in an earlier surface "
            "table, returning only the current coupling-window increment"
        ),
    )
    return_parser.add_argument(
        "--polar-wander-history",
        type=Path,
        help=(
            "map body-fixed topography into the spin frame defined by the "
            "last pole in this true_polar_wander.csv file"
        ),
    )
    return_parser.add_argument(
        "--interpolation-neighbors",
        type=int,
        default=4,
        help=(
            "number of nearby Fastscape cells used for continuous spherical "
            "weighting; one recovers the original nearest-cell mapping"
        ),
    )
    arguments = parser.parse_args()

    if arguments.command == "inspect":
        exchange = read_exchange(arguments.exchange)
        print(f"points={exchange.longitude.size} time_years={exchange.model_time_years:g}")
        for name, values in exchange.fields.items():
            print(f"{name} [{exchange.units[name]}]: min={values.min():.9g} max={values.max():.9g}")
    elif arguments.command == "climate-to-aspect":
        exchange = read_exchange(arguments.exchange)
        erosion_strength, runoff = climate_controls(exchange)
        if arguments.surface_topography is not None:
            write_aspect_structured(
                arguments.surface_topography,
                exchange,
                exchange.fields["surface_elevation"],
            )
        write_aspect_structured(arguments.erosion_strength, exchange, erosion_strength)
        write_aspect_structured(arguments.surface_runoff, exchange, runoff)
        if arguments.ice_thickness is not None:
            write_aspect_structured(
                arguments.ice_thickness, exchange, exchange.fields["ice_thickness"]
            )
        if arguments.basal_ice_velocity is not None:
            basal_ice_velocity = exchange.fields["basal_ice_velocity"]
            if arguments.prescribed_basal_ice_velocity is not None:
                if arguments.prescribed_basal_ice_velocity < 0.0:
                    parser.error("prescribed basal ice velocity must be nonnegative")
                basal_ice_velocity = np.where(
                    exchange.fields["grounded_ice_fraction"] > 0.0,
                    arguments.prescribed_basal_ice_velocity,
                    0.0,
                )
            write_aspect_structured(
                arguments.basal_ice_velocity,
                exchange,
                basal_ice_velocity,
            )
    else:
        spin_axis = (
            spin_axis_from_history(arguments.polar_wander_history)
            if arguments.polar_wander_history is not None
            else None
        )
        surface_to_climate_exchange(
            arguments.surface,
            arguments.climate,
            arguments.output,
            spin_axis,
            arguments.interpolation_neighbors,
            arguments.reference_surface,
        )


if __name__ == "__main__":
    main()
