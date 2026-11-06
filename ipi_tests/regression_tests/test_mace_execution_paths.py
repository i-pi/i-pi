"""End-to-end regression test for the three MACE execution paths."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from ase.io import read

pytest.importorskip("mace")

REPO_ROOT = Path(__file__).parents[2]
EXAMPLE_ROOT = REPO_ROOT / "examples" / "clients" / "mace" / "single"
COMPARISON_ROOT = EXAMPLE_ROOT / "comparison"
POSITION_ATOL = 1e-9
PROPERTY_ATOL = 1e-9


@pytest.fixture
def mace_mp_model_path():
    """Return the optional model downloaded by the MACE example."""

    model_path = EXAMPLE_ROOT / "mace.model"
    if not model_path.is_file():
        pytest.skip(
            "requires the MACE-MP model downloaded by "
            "examples/clients/mace/single/getmodel.sh"
        )
    return model_path


def _positions(path):
    return np.asarray([atoms.positions for atoms in read(path, index=":")])


def _failure_details(results):
    paths = [results / "FAILED", results / "run.log"]
    paths.extend(sorted(results.glob("*.log")))
    details = []
    for path in paths:
        if path.is_file():
            details.append(f"--- {path.name}\n{path.read_text(errors='replace')}")
    return "\n".join(details)


def test_mace_execution_paths_produce_same_dynamics(tmp_path, mace_mp_model_path):
    """Run socket clients and ffdirect, then compare their NVE trajectories."""

    workdir = tmp_path / "single" / "comparison"
    shutil.copytree(
        COMPARISON_ROOT,
        workdir,
        ignore=shutil.ignore_patterns(
            "__pycache__", "results", "results.previous.*", "error.txt"
        ),
    )
    (workdir.parent / "mace.model").symlink_to(mace_mp_model_path)

    environment = os.environ.copy()
    environment["PYTHON"] = sys.executable
    environment["IPI_REPO_ROOT"] = str(REPO_ROOT)
    completed = subprocess.run(
        ["bash", "run.sh"],
        cwd=workdir,
        env=environment,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )

    results = workdir / "results"
    if completed.returncode != 0:
        pytest.fail(
            "MACE execution-path comparison failed:\n"
            f"{completed.stdout}\n{completed.stderr}\n"
            f"{_failure_details(results)}"
        )

    assert (results / "SUCCESS").is_file()

    trajectory_names = {
        "macecalculator": "ffsocket_macecalculator.trajectory_0.xyz",
        "py_driver": "ffsocket_py_driver.trajectory_0.xyz",
        "ffdirect": "ffdirect.trajectory_0.xyz",
    }
    property_names = {
        key: filename.replace("trajectory_0.xyz", "properties")
        for key, filename in trajectory_names.items()
    }
    reference_positions = _positions(results / trajectory_names["ffdirect"])
    reference_properties = np.loadtxt(results / property_names["ffdirect"])

    for key in ("macecalculator", "py_driver"):
        np.testing.assert_allclose(
            _positions(results / trajectory_names[key]),
            reference_positions,
            rtol=1e-9,
            atol=POSITION_ATOL,
        )
        np.testing.assert_allclose(
            np.loadtxt(results / property_names[key]),
            reference_properties,
            rtol=1e-9,
            atol=PROPERTY_ATOL,
        )
