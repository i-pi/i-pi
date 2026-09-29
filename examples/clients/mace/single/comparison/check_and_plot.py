#!/usr/bin/env python3
"""Check three MACE trajectories and plot their end-to-end wall times."""

import csv
from pathlib import Path

import matplotlib
import numpy as np
from ase.io import read

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = Path("results")
RELATIVE_TOLERANCE = 1e-9
ABSOLUTE_TOLERANCE = 1e-9
CASES = {
    "ffsocket_macecalculator": {
        "label": "ffsocket +\nMACECalculator",
        "trajectory": RESULTS / "ffsocket_macecalculator.trajectory_0.xyz",
        "properties": RESULTS / "ffsocket_macecalculator.properties",
    },
    "ffsocket_py_driver": {
        "label": "ffsocket +\ni-pi-py_driver",
        "trajectory": RESULTS / "ffsocket_py_driver.trajectory_0.xyz",
        "properties": RESULTS / "ffsocket_py_driver.properties",
    },
    "ffdirect": {
        "label": "ffdirect",
        "trajectory": RESULTS / "ffdirect.trajectory_0.xyz",
        "properties": RESULTS / "ffdirect.properties",
    },
}


def load_trajectory(path):
    if not path.is_file():
        raise FileNotFoundError(f"missing trajectory: {path}")
    return read(path, index=":")


def load_properties(path):
    if not path.is_file():
        raise FileNotFoundError(f"missing properties file: {path}")
    values = np.loadtxt(path)
    return np.atleast_2d(values)


def check_case(
    reference,
    candidate,
    name,
    rtol=RELATIVE_TOLERANCE,
    atol=ABSOLUTE_TOLERANCE,
):
    reference_frames = load_trajectory(reference["trajectory"])
    candidate_frames = load_trajectory(candidate["trajectory"])
    if len(reference_frames) != len(candidate_frames):
        raise AssertionError(
            f"{name}: found {len(candidate_frames)} frames; "
            f"expected {len(reference_frames)}"
        )

    reference_positions = np.asarray([frame.positions for frame in reference_frames])
    candidate_positions = np.asarray([frame.positions for frame in candidate_frames])
    for index, (reference_frame, candidate_frame) in enumerate(
        zip(reference_frames, candidate_frames)
    ):
        if not np.array_equal(reference_frame.numbers, candidate_frame.numbers):
            raise AssertionError(f"{name}: atomic species differ in frame {index}")
        np.testing.assert_allclose(
            reference_frame.cell.array,
            candidate_frame.cell.array,
            rtol=rtol,
            atol=atol,
            err_msg=f"{name}: cell differs in frame {index}",
        )

    np.testing.assert_allclose(
        reference_positions,
        candidate_positions,
        rtol=rtol,
        atol=atol,
        err_msg=f"{name}: coordinates differ",
    )

    reference_properties = load_properties(reference["properties"])
    candidate_properties = load_properties(candidate["properties"])
    np.testing.assert_allclose(
        reference_properties,
        candidate_properties,
        rtol=rtol,
        atol=atol,
        err_msg=f"{name}: energies or conserved quantities differ",
    )

    coordinate_error = np.max(np.abs(reference_positions - candidate_positions))
    property_error = np.max(np.abs(reference_properties - candidate_properties))
    print(
        f"PASS: {name}: {len(reference_frames) - 1} NVE steps; "
        f"max |Δposition| = {coordinate_error:.3e} Å; "
        f"max |Δproperty| = {property_error:.3e}"
    )


def load_timings():
    with open(RESULTS / "timings.csv", newline="", encoding="utf-8") as handle:
        timings = {row["case"]: float(row["seconds"]) for row in csv.DictReader(handle)}
    missing = set(CASES) - set(timings)
    if missing:
        raise ValueError(f"missing timings for: {', '.join(sorted(missing))}")
    return timings


def make_plot(timings):
    keys = list(CASES)
    labels = [CASES[key]["label"] for key in keys]
    seconds = [timings[key] for key in keys]
    fastest = int(np.argmin(seconds))
    colors = ["#4C78A8", "#F58518", "#54A24B"]
    colors[fastest] = "#2CA02C"

    figure, axis = plt.subplots(figsize=(8, 4.8), layout="constrained")
    bars = axis.bar(labels, seconds, color=colors)
    axis.set_ylabel("End-to-end wall time (s)")
    axis.set_title("MACE execution-path timing (10 NVE steps)")
    axis.grid(axis="y", alpha=0.25)
    axis.set_axisbelow(True)
    axis.bar_label(bars, labels=[f"{value:.2f} s" for value in seconds], padding=3)
    axis.text(
        0.99,
        0.97,
        f"Fastest: {CASES[keys[fastest]]['label'].replace(chr(10), ' ')}",
        transform=axis.transAxes,
        ha="right",
        va="top",
    )
    figure.savefig(RESULTS / "timing.png", dpi=180)
    plt.close(figure)
    print(f"PASS: wrote {RESULTS / 'timing.png'}")


def main():
    reference = CASES["ffdirect"]
    check_case(
        reference,
        CASES["ffsocket_macecalculator"],
        "ffsocket + MACECalculator versus ffdirect",
    )
    check_case(
        reference,
        CASES["ffsocket_py_driver"],
        "ffsocket + i-pi-py_driver versus ffdirect",
    )
    make_plot(load_timings())


if __name__ == "__main__":
    main()
