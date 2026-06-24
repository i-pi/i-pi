"""End-to-end FFMPI test for the q-TIP4P/f water example.

Complements test_ffmpi_mpirun.py (harmonic) and test_ffmpi_multi_mpirun.py
(two tagged drivers) by covering the remaining communication/mpi case: the compiled
Fortran driver returning forces AND a dipole `extras` string over MPI, with the
ring-polymer beads sent as one batched request. Needs mpi4py, a working MPI
runtime (IPI_TEST_MPI=1) and the Fortran driver built with `make MPI=1`.
"""

import os
import shutil
import subprocess
import importlib.util
from pathlib import Path

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
    / "water-dipole"
)
DRIVER = Path(__file__).resolve().parents[3] / "drivers" / "f90" / "driver.x"


def _fortran_mpi_built():
    """True if the compiled Fortran driver links against MPI."""
    if not DRIVER.exists():
        return False
    try:
        out = subprocess.run(["nm", str(DRIVER)], capture_output=True, text=True).stdout
    except OSError:
        return False
    return "mpi_init" in out.lower()


@pytest.mark.skipif(
    not _fortran_mpi_built(),
    reason="needs the q-TIP4P/f Fortran driver built with `make MPI=1`",
)
def test_ffmpi_water_mpirun(tmp_path):
    mpirun = shutil.which("mpirun")
    assert mpirun is not None, "mpirun not found on PATH"

    for name in ("input.xml", "water.xyz"):
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
        "1",
        str(DRIVER),
        "--mpi",
        "-m",
        "qtip4pf",
    ]
    proc = subprocess.run(
        cmd, cwd=tmp_path, capture_output=True, text=True, timeout=300
    )
    assert proc.returncode == 0, proc.stderr

    # the dipole is an 'extras' field returned by the driver over MPI; its
    # trajectory file is only written if the extras round-trip works
    assert (tmp_path / "simulation.dipole_0").exists(), "dipole extras were not written"
