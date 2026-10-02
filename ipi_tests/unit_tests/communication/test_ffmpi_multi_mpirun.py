"""End-to-end test for MULTIPLE ffmpi forcefields in one MPMD job.

Needs mpi4py + a working MPI runtime, so skipped unless mpi4py is importable and
IPI_TEST_MPI=1. Launches i-PI with two `<ffmpi>` forcefields, each served by its
own driver tagged with --mpi-name, and checks the trajectory is bit-identical to
the equivalent run with two in-process `<ffdirect>` harmonic forcefields (which
confirms the driver ranks were routed to the right forcefield).
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
    / "two-drivers"
)

# the two ffdirect blocks equivalent to the tagged harmonic drivers (k=1 and k=2)
_FFDIRECT = (
    "<ffdirect name='ffA'><pes>harmonic</pes>"
    "<parameters>{k1:1.0}</parameters></ffdirect>\n"
    "<ffdirect name='ffB'><pes>harmonic</pes>"
    "<parameters>{k1:2.0}</parameters></ffdirect>"
)


def test_ffmpi_multi_mpirun(tmp_path):
    mpirun = shutil.which("mpirun")
    assert mpirun is not None, "mpirun not found on PATH"

    shutil.copy(EXAMPLE / "init.xyz", tmp_path / "init.xyz")
    shutil.copy(EXAMPLE / "input.xml", tmp_path / "input.xml")

    # reference: replace the two <ffmpi> blocks with equivalent <ffdirect> ones
    ref_xml = (tmp_path / "input.xml").read_text()
    ref_xml = (
        ref_xml.replace("<ffmpi name='ffA' mode='mpi'></ffmpi>", "")
        .replace("<ffmpi name='ffB' mode='mpi'></ffmpi>", _FFDIRECT)
        .replace("prefix='simulation'", "prefix='ref'")
    )
    (tmp_path / "ref.xml").write_text(ref_xml)

    ref = subprocess.run(
        ["i-pi", "ref.xml"], cwd=tmp_path, capture_output=True, text=True, timeout=300
    )
    assert ref.returncode == 0, ref.stderr

    cmd = [mpirun, "--oversubscribe", "-n", "1", "i-pi", "input.xml"]
    for name, k in (("ffA", "1"), ("ffB", "2")):
        cmd += [
            ":",
            "-n",
            "1",
            "i-pi-py_driver",
            "--mpi",
            "--mpi-name",
            name,
            "-m",
            "harmonic",
            "-o",
            k,
        ]
    multi = subprocess.run(
        cmd, cwd=tmp_path, capture_output=True, text=True, timeout=300
    )
    assert multi.returncode == 0, multi.stderr

    got = (tmp_path / "simulation.pos_0.xyz").read_text()
    want = (tmp_path / "ref.pos_0.xyz").read_text()
    assert got == want, "two-ffmpi trajectory differs from the two-ffdirect reference"
