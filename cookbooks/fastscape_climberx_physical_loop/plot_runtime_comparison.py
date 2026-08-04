#!/usr/bin/env python3
"""Plot measured coupled and uncoupled runtimes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    data = json.loads(arguments.summary.read_text(encoding="utf-8"))
    coupled_aspect = data["coupled_aspect_seconds"]
    climate = data["climate_seconds"]
    exchange = sum(data["exchange_seconds"])
    overhead = max(
        data["full_coupling_seconds"] - coupled_aspect - climate - exchange, 0
    )
    uncoupled = data["uncoupled_aspect_seconds"]

    figure, axis = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    bottom = 0.0
    for value, label, color in (
        (coupled_aspect, "ASPECT and Fastscape", "#3973ac"),
        (climate, "CLIMBER-X and Yelmo", "#df8f2d"),
        (exchange, "field interpolation", "#4ca36b"),
        (overhead, "file and process overhead", "#999999"),
    ):
        axis.bar(0, value, bottom=bottom, label=label, color=color)
        bottom += value
    axis.bar(1, uncoupled, color="#3973ac")
    ticks = [0, 1]
    labels = ["full coupling", "ASPECT only\ncontinuous"]
    if "uncoupled_windowed_total_seconds" in data:
        axis.bar(2, data["uncoupled_windowed_total_seconds"], color="#3973ac")
        ticks.append(2)
        labels.append("ASPECT only\nthree windows")
    axis.set_xticks(ticks, labels)
    axis.set_ylabel("measured wall time (seconds)")
    axis.set_title(
        f"{data['coupling_windows']} coupling windows; "
        f"full loop is {data['full_to_uncoupled_ratio']:.1f}× slower"
    )
    axis.legend(frameon=False, loc="upper right")
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output, dpi=180)


if __name__ == "__main__":
    main()
