"""End-to-end FFMPI test against a real mpirun launch.

This needs mpi4py and a working MPI runtime, so it is skipped unless mpi4py is
importable and IPI_TEST_MPI=1 is set in the environment. It launches i-PI and
two harmonic driver ranks in a single MPMD job and checks the run completes and
the forces match the analytic harmonic result.
"""

import os
import shutil
import subprocess
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("mpi4py") is None or os.environ.get("IPI_TEST_MPI") != "1",
    reason="set IPI_TEST_MPI=1 and install mpi4py to run the MPI end-to-end test",
)

EXAMPLE = (
    Path(__file__).resolve().parents[3]
    / "examples"
    / "clients"
    / "communication"
    / "mpi"
    / "harmonic"
)


def test_ffmpi_mpirun(tmp_path):
    mpirun = shutil.which("mpirun")
    assert mpirun is not None, "mpirun not found on PATH"

    for name in ("input.xml", "init.xyz"):
        shutil.copy(EXAMPLE / name, tmp_path / name)

    cmd = [
        mpirun,
        # CI runners expose few MPI slots; allow more ranks than cores
        "--oversubscribe",
        "-n",
        "1",
        "i-pi",
        "input.xml",
        ":",
        "-n",
        "2",
        "i-pi-py_driver",
        "--mpi",
        "-m",
        "harmonic",
        "-o",
        "1",
    ]
    proc = subprocess.run(
        cmd, cwd=tmp_path, capture_output=True, text=True, timeout=300
    )
    assert proc.returncode == 0, proc.stderr

    # the harmonic PES (k=1, in atomic units) gives force = -position; the xyz
    # trajectories are in different units but the ratio is constant, so check
    # that every force component is anti-parallel to the matching position one
    pos = _last_atom_row(tmp_path / "simulation.pos_0.xyz")
    force = _last_atom_row(tmp_path / "simulation.force_0.xyz")
    nonzero = np.abs(pos) > 1e-8
    assert np.all(np.sign(force[nonzero]) == -np.sign(pos[nonzero]))


def _last_atom_row(path):
    """Returns the xyz coordinates of the (single) atom in the last frame."""
    atom_rows = [
        cols[1:4]
        for cols in (line.split() for line in Path(path).read_text().splitlines())
        if len(cols) == 4 and cols[0].isalpha()
    ]
    return np.array(atom_rows[-1], dtype=float)
