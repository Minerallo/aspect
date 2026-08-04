#!/usr/bin/env python3
"""Continue several coupling windows and benchmark an uncoupled ASPECT run."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

from run_physical_loop import (
    ASPECT_SOURCE,
    COOKBOOK,
    EXCHANGE_TOOL,
    INSTALLATION,
    convert_climate,
    newest,
    prepare_climate_case,
    run_climate,
    run_monitored,
    write_aspect_parameters,
)


def replace_namelist_value(text: str, name: str, value: str) -> str:
    pattern = re.compile(rf"^(\s*{re.escape(name)}\s*=\s*)[^!\n]*(.*)$", re.MULTILINE)
    replaced, count = pattern.subn(rf"\g<1>{value} \g<2>", text, count=1)
    if count != 1:
        raise RuntimeError(f"could not set {name} in control.nml")
    return replaced


def configure_continuation(case: Path, restart_directory: Path, year: int) -> None:
    control = case / "control.nml"
    text = control.read_text(encoding="utf-8")
    text = replace_namelist_value(text, "year_ini", str(year))
    text = replace_namelist_value(text, "nyears", "1")
    for name in (
        "atm_restart",
        "lnd_restart",
        "ocn_restart",
        "sic_restart",
        "geo_restart",
        "ice_restart",
        "smb_restart",
        "bmb_restart",
    ):
        text = replace_namelist_value(text, name, "true")
    text = replace_namelist_value(
        text, "restart_in_dir", f'"{restart_directory.resolve()}"'
    )
    control.write_text(text, encoding="utf-8")


def write_uncoupled_parameters(
    path: Path,
    output_directory: Path,
    topography_directory: Path,
    end_time: float,
    resume: str = "false",
) -> None:
    base = (
        INSTALLATION
        / "aspect_fatscapecc/tests/fastscape_cpp_global_climate_geodynamics.prm"
    ).resolve()
    path.write_text(
        f"""# Matched ASPECT benchmark without the surface-process coupling.
include {base}

set Output directory   = {output_directory.resolve()}
set Resume computation = {resume}
set End time            = {end_time:g}
set Maximum time step   = 10

subsection Geometry model
  subsection Initial topography model
    set Model name = ascii data
    subsection Ascii data model
      set Data directory = {topography_directory.resolve()}/
      set Data file name = surface-topography.txt
    end
  end
end

subsection Mesh deformation
  set Mesh deformation boundary indicators =
end
""",
        encoding="utf-8",
    )


def return_surface_increment(
    surface: Path,
    climate_exchange: Path,
    output: Path,
    polar_history: Path | None,
    reference_surface: Path | None,
) -> float:
    command = [
        sys.executable,
        str(EXCHANGE_TOOL),
        "surface-to-climate",
        "--surface",
        str(surface),
        "--climate",
        str(climate_exchange),
        "--output",
        str(output),
    ]
    if polar_history is not None:
        command.extend(("--polar-wander-history", str(polar_history)))
    if reference_surface is not None:
        command.extend(("--reference-surface", str(reference_surface)))
    start = time.perf_counter()
    subprocess.run(command, check=True)
    return time.perf_counter() - start


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", type=int, default=3)
    parser.add_argument(
        "--output", type=Path, default=COOKBOOK / "output-coupling-sequence"
    )
    parser.add_argument("--minimum-free-memory-percent", type=int, default=15)
    parser.add_argument(
        "--windowed-benchmark-only",
        action="store_true",
        help="add the matched checkpointed ASPECT-only benchmark to existing results",
    )
    arguments = parser.parse_args()
    if arguments.windows < 2:
        parser.error("use at least two coupling windows for a sequence")

    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    climate_root = INSTALLATION / "aspect_ClimberX"
    climate_template = climate_root / "tests/yelmo-active-diva"
    climate_executable = climate_root / "sources/climber-x/climber.x"
    aspect_executable = (
        INSTALLATION / "aspect_fatscapecc/builts/fastscape-release/aspect-release"
    )
    if arguments.windowed_benchmark_only:
        summary_path = output / "timing-summary.json"
        if not summary_path.is_file():
            parser.error(f"existing timing summary is missing: {summary_path}")
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        windowed_output = output / "aspect-uncoupled-windowed"
        windowed_seconds = []
        for window in range(1, arguments.windows + 1):
            parameter_file = output / f"aspect-uncoupled-window-{window:03d}.prm"
            write_uncoupled_parameters(
                parameter_file,
                windowed_output,
                output / "aspect-input-001",
                20.0 * window,
                resume="false" if window == 1 else "auto",
            )
            windowed_seconds.append(
                run_monitored(
                    [str(aspect_executable), str(parameter_file)],
                    INSTALLATION / "aspect_fatscapecc/tests",
                    output / f"aspect-uncoupled-window-{window:03d}.log",
                    os.environ.copy(),
                    arguments.minimum_free_memory_percent,
                )
            )
        summary["uncoupled_windowed_seconds"] = windowed_seconds
        summary["uncoupled_windowed_total_seconds"] = sum(windowed_seconds)
        summary["coupled_aspect_to_windowed_uncoupled_ratio"] = summary[
            "coupled_aspect_seconds"
        ] / sum(windowed_seconds)
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        plot_environment = os.environ.copy()
        plot_environment["MPLCONFIGDIR"] = str(output / ".matplotlib")
        subprocess.run(
            [
                sys.executable,
                str(COOKBOOK / "plot_runtime_comparison.py"),
                str(summary_path),
                "--output",
                str(output / "runtime-comparison.png"),
            ],
            check=True,
            env=plot_environment,
        )
        print(json.dumps(summary, indent=2))
        return
    total_start = time.perf_counter()
    climate_seconds = []
    aspect_seconds = []
    exchange_seconds = []

    climate_case = output / "climate-000"
    prepare_climate_case(climate_template, climate_case, climate_executable)
    climate_exchange = output / "climate-000.cxe"
    climate_seconds.append(
        run_climate(
            climate_case,
            climate_exchange,
            None,
            climate_root / "local/lib",
            arguments.minimum_free_memory_percent,
        )
    )

    aspect_output = output / "aspect-fastscape"
    previous_surface = None
    polar_history = None
    for window in range(1, arguments.windows + 1):
        fields = output / f"aspect-input-{window:03d}"
        conversion_start = time.perf_counter()
        convert_climate(climate_exchange, fields)
        exchange_seconds.append(time.perf_counter() - conversion_start)
        parameter_file = output / f"aspect-window-{window:03d}.prm"
        write_aspect_parameters(
            parameter_file,
            aspect_output,
            fields,
            end_time=20.0 * window,
            resume="false" if window == 1 else "auto",
        )
        aspect_seconds.append(
            run_monitored(
                [str(aspect_executable), str(parameter_file)],
                INSTALLATION / "aspect_fatscapecc/tests",
                output / f"aspect-window-{window:03d}.log",
                os.environ.copy(),
                arguments.minimum_free_memory_percent,
            )
        )
        current_surface = newest(aspect_output, "surface-*.csv")
        polar_history = next(iter(aspect_output.rglob("true_polar_wander.csv")), None)
        topography_exchange = output / f"topography-{window:03d}.cxe"
        exchange_seconds.append(
            return_surface_increment(
                current_surface,
                climate_exchange,
                topography_exchange,
                polar_history,
                previous_surface,
            )
        )

        next_case = output / f"climate-{window:03d}"
        prepare_climate_case(climate_template, next_case, climate_executable)
        restart_directory = climate_case / "restart_out" / f"year_{window}"
        if not restart_directory.is_dir():
            raise RuntimeError(f"missing climate restart: {restart_directory}")
        configure_continuation(next_case, restart_directory, window)
        next_exchange = output / f"climate-{window:03d}.cxe"
        climate_seconds.append(
            run_climate(
                next_case,
                next_exchange,
                topography_exchange,
                climate_root / "local/lib",
                arguments.minimum_free_memory_percent,
            )
        )
        previous_surface = current_surface
        climate_case = next_case
        climate_exchange = next_exchange

    coupled_total_seconds = time.perf_counter() - total_start
    uncoupled_parameters = output / "aspect-uncoupled.prm"
    uncoupled_output = output / "aspect-uncoupled"
    write_uncoupled_parameters(
        uncoupled_parameters,
        uncoupled_output,
        output / "aspect-input-001",
        20.0 * arguments.windows,
    )
    uncoupled_seconds = run_monitored(
        [str(aspect_executable), str(uncoupled_parameters)],
        INSTALLATION / "aspect_fatscapecc/tests",
        output / "aspect-uncoupled.log",
        os.environ.copy(),
        arguments.minimum_free_memory_percent,
    )

    summary = {
        "coupling_windows": arguments.windows,
        "surface_model_years": 20 * arguments.windows,
        "climate_model_years": arguments.windows + 1,
        "climate_window_seconds": climate_seconds,
        "coupled_aspect_window_seconds": aspect_seconds,
        "exchange_seconds": exchange_seconds,
        "full_coupling_seconds": coupled_total_seconds,
        "coupled_aspect_seconds": sum(aspect_seconds),
        "climate_seconds": sum(climate_seconds),
        "uncoupled_aspect_seconds": uncoupled_seconds,
        "full_to_uncoupled_ratio": coupled_total_seconds / uncoupled_seconds,
        "coupled_aspect_to_uncoupled_ratio": sum(aspect_seconds) / uncoupled_seconds,
        "final_climate_exchange": str(climate_exchange),
        "final_surface": str(previous_surface),
        "polar_wander_history": str(polar_history) if polar_history else None,
    }
    summary_path = output / "timing-summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    plot_environment = os.environ.copy()
    plot_environment["MPLCONFIGDIR"] = str(output / ".matplotlib")
    subprocess.run(
        [
            sys.executable,
            str(COOKBOOK / "plot_runtime_comparison.py"),
            str(summary_path),
            "--output",
            str(output / "runtime-comparison.png"),
        ],
        check=True,
        env=plot_environment,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
