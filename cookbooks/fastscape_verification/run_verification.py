#!/usr/bin/env python3
"""Run the bounded-memory FastScape and CLIMBER-X verification workflow."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import psutil


HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[3]


def run_monitored(command: list[str], cwd: Path, log: Path,
                  minimum_free_fraction: float = 0.15) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with log.open("w", encoding="utf-8") as stream:
        process = subprocess.Popen(command, cwd=cwd, stdout=stream,
                                   stderr=subprocess.STDOUT)
        while process.poll() is None:
            memory = psutil.virtual_memory()
            if memory.available / memory.total < minimum_free_fraction:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                raise RuntimeError(
                    f"Stopped {' '.join(command)} because free memory fell below "
                    f"{100 * minimum_free_fraction:.0f}%"
                )
            time.sleep(1)
    if process.returncode:
        raise RuntimeError(f"Command failed; see {log}")
    print(f"PASS {command[-1]} ({time.monotonic() - started:.2f} s)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aspect-executable", type=Path, default=WORKSPACE /
                        "aspect_fatscapecc/builts/fastscape-release/aspect-release")
    parser.add_argument("--fortran-root", type=Path, default=WORKSPACE /
                        "aspect_fastscapef90_glacial")
    parser.add_argument("--cpp-build", type=Path, default=WORKSPACE /
                        "aspect_fatscapecc/builts/test_newfeatures")
    parser.add_argument("--spherical-build", type=Path, default=WORKSPACE /
                        "aspect_sphericaldiffusion/builts")
    parser.add_argument("--external-tests", type=Path, default=WORKSPACE /
                        "aspect_fatscapecc/tests")
    parser.add_argument("--analysis-only", action="store_true")
    parser.add_argument("--run-two-process", action="store_true",
                        help="Repeat the local two-process message-passing tests.")
    parser.add_argument("--run-climate", action="store_true",
                        help="Repeat the several-minute CLIMBER-X physical loop.")
    args = parser.parse_args()

    output = HERE / "output"
    logs = output / "logs"
    output.mkdir(exist_ok=True)
    os.environ.setdefault("MACOSX_DEPLOYMENT_TARGET", "13.0")
    os.environ.setdefault("DEVELOPER_DIR", "/Library/Developer/CommandLineTools")
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/fastscape-matplotlib")
    environment = os.environ.copy()

    if not args.analysis_only:
        compiler = shutil.which("mpif90") or shutil.which("gfortran")
        if compiler is None:
            raise RuntimeError("No Fortran compiler was found")
        library = args.fortran_root / "build-glacial-release/libfastscapelib_fortran.a"
        executable = HERE / "build/fortran_reference_benchmarks"
        executable.parent.mkdir(exist_ok=True)
        compile_command = [
            compiler, "-cpp", "-ffree-form", "-std=f2008", "-fimplicit-none",
            "-fall-intrinsics", "-fconvert=big-endian", "-O2",
            str(HERE / "fortran_reference_benchmarks.f90"), "-o", str(executable),
            str(library),
        ]
        subprocess.run(compile_command, cwd=HERE, env=environment, check=True)
        run_monitored([str(executable)], HERE, logs / "fortran-reference.log")

        run_monitored(["ctest", "--test-dir", str(args.cpp_build),
                       "--output-on-failure", "-j2"], HERE,
                      logs / "fastscapelib-cpp-tests.log")
        run_monitored(["ctest", "--test-dir",
                       str(args.fortran_root / "build-glacial-release"),
                       "--output-on-failure", "-j2"], HERE,
                      logs / "fastscape-fortran-tests.log")
        run_monitored(["ctest", "--test-dir", str(args.spherical_build),
                       "--output-on-failure", "-j2", "-R",
                       "^hill_diffusion_spherical_(2d|3d)$"], HERE,
                      logs / "spherical-diffusion-tests.log")

        unit_tests = args.aspect_executable.parent / "bin/unit_tests"
        run_monitored([str(unit_tests)], HERE, logs / "aspect-unit-tests.log")
        exchange_test = (
            HERE.parents[1] / "contrib/coupling/climberx_fastscape"
            / "test_surface_exchange.py"
        )
        run_monitored([sys.executable, "-m", "pytest", "-q", str(exchange_test)],
                      HERE, logs / "climber-x-exchange-tests.log")

        serial_parameters = (
            "fastscape_cpp_box.prm",
            "fastscape_cpp_box_glacial_erosion.prm",
            "fastscape_cpp_global.prm",
            "fastscape_cpp_global_strong.prm",
            "fastscape_cpp_global_zero_topography.prm",
            "fastscape_cpp_global_gmg.prm",
            "fastscape_cpp_global_marine_deposition.prm",
            "fastscape_cpp_global_advection.prm",
            "fastscape_cpp_global_climber_glacial.prm",
        )
        for parameter_file in serial_parameters:
            run_monitored([str(args.aspect_executable), parameter_file],
                          args.external_tests,
                          logs / f"{parameter_file}.log")

        if args.run_two_process:
            message_passing_runner = shutil.which("mpirun")
            if message_passing_runner is None:
                raise RuntimeError("The message-passing launcher was not found")
            for parameter_file in (
                "fastscape_cpp_box_mpi.prm",
                "fastscape_cpp_global_mpi.prm",
            ):
                run_monitored([message_passing_runner, "-np", "2",
                               str(args.aspect_executable), parameter_file],
                              args.external_tests,
                              logs / f"{parameter_file}.log")

        for validator in (
            "validate_outputs.py",
            "validate_global_suite.py",
            "validate_marine_deposition.py",
        ):
            run_monitored([sys.executable, validator], args.external_tests,
                          logs / f"{validator}.log")

        for parameter_file in (
            "global_flat_surface_control.prm",
            "global_rigid_rotation_control.prm",
            "global_rigid_rotation.prm",
            "global_rigid_rotation_refined.prm",
            "global_rigid_rotation_fine.prm",
        ):
            run_monitored([str(args.aspect_executable), parameter_file], HERE,
                          logs / f"{parameter_file}.log")

        if args.run_climate:
            climate_runner = (
                HERE.parent / "fastscape_climberx_physical_loop"
                / "run_coupling_sequence.py"
            )
            run_monitored([sys.executable, str(climate_runner),
                           "--ice-model", "diagnostic"], climate_runner.parent,
                          logs / "climber-x-physical-loop.log")

    run_monitored(
        [sys.executable, str(HERE / "analyze_verification.py"),
         "--external-tests", str(args.external_tests)],
        HERE, logs / "analysis.log"
    )


if __name__ == "__main__":
    main()
