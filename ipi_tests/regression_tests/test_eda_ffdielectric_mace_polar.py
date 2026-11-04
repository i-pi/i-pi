"""Compare legacy EDA dynamics with the FFDielectric implementation."""

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ase")
pytest.importorskip("mace")

REPO_ROOT = Path(__file__).parents[2]
MODEL_CANDIDATES = (
    REPO_ROOT
    / "examples/features/ffdieletric/mace-polar+E-server/MACE-POLAR-1-M.model",
    REPO_ROOT / "examples/clients/mace/polar/MACE-POLAR-1-M.model",
)
SOURCE_EXAMPLE = REPO_ROOT / "examples/features/ffdieletric/mace-polar+E-server"
BEC_SOURCE = REPO_ROOT / "ipi_tests/regression_tests/tests/EDANVE/bec.txt"
DRIVER = Path(__file__).with_name("mace_polar_response_driver.py")


@pytest.fixture
def mace_polar_model():
    for candidate in MODEL_CANDIDATES:
        if candidate.is_file():
            return candidate
    pytest.skip(
        "MACE-POLAR model not found; run "
        "examples/features/ffdieletric/mace-polar+E-server/download_mace_polar.sh"
    )


def _wait_for_socket(address, process):
    socket_path = Path("/tmp") / f"ipi_{address}"
    for _ in range(300):
        if socket_path.exists():
            return
        if process.poll() is not None:
            raise RuntimeError("i-PI exited before creating the socket")
        time.sleep(0.1)
    raise RuntimeError(f"Timed out waiting for {socket_path}")


def _field_file(path):
    path.write_text(
        "import numpy as np\n\n"
        "def field(time, amplitude, frequency, peak, sigma):\n"
        "    envelope = np.exp(-0.5 * ((time - peak) / sigma) ** 2)\n"
        "    return np.asarray(amplitude) * np.cos(frequency * time) * envelope\n"
    )


def _input_xml(implementation, ensemble, address):
    if implementation == "eda":
        forcefield = f"""
  <ffsocket name='driver' mode='unix' pbc='false'>
    <address>{address}</address>
  </ffsocket>"""
        motion = f"""
    <motion mode='driven_dynamics'>
      <driven_dynamics mode='eda-{ensemble}'>
        <timestep units='femtosecond'>0.5</timestep>
        <efield>
          <amp units='atomic_unit'>[0.0, 0.002, 0.0]</amp>
          <freq units='atomic_unit'>0.001</freq>
          <peak units='atomic_unit'>2.0</peak>
          <sigma units='atomic_unit'>4.0</sigma>
        </efield>
        <bec mode='file'>{BEC_SOURCE.name}</bec>
      </driven_dynamics>
    </motion>"""
    else:
        forcefield = f"""
  <ffdielectric name='mu+E'>
    <electric_field file='field.py' name='field' units='atomic_unit'>
      <parameters>{{"amplitude": [0.0, 0.002, 0.0], "frequency": 0.001,
                    "peak": 2.0, "sigma": 4.0}}</parameters>
    </electric_field>
    <ffsocket name='driver' mode='unix' pbc='false'>
      <address>{address}</address>
    </ffsocket>
  </ffdielectric>"""
        motion = f"""
    <motion mode='dynamics'>
      <dynamics mode='{ensemble}'>
        <timestep units='femtosecond'>0.5</timestep>
        {"<thermostat mode='langevin'><tau units='femtosecond'>10.0</tau></thermostat>" if ensemble == "nvt" else ""}
      </dynamics>
    </motion>"""

    return f"""<simulation verbosity='low'>
  <total_steps>6</total_steps>
  <prng><seed>31415</seed></prng>
{forcefield}
  <output prefix='simulation'>
    <trajectory filename='positions' format='ase'>positions</trajectory>
    <trajectory filename='velocities' format='ase'>velocities</trajectory>
  </output>
  <system>
    <initialize nbeads='1'>
      <file mode='ase'>start.extxyz</file>
      <velocities mode='manual' units='atomic_unit'>[0, 0, 0, 0, 0, 0, 0, 0, 0]</velocities>
    </initialize>
    <forces><force forcefield='{'driver' if implementation == 'eda' else 'mu+E'}'/></forces>
    <ensemble><temperature units='kelvin'>0</temperature></ensemble>
{motion}
  </system>
</simulation>
"""


def _run_case(tmp_path, implementation, ensemble, model):
    case = tmp_path / f"{implementation}-{ensemble}"
    case.mkdir()
    shutil.copy(SOURCE_EXAMPLE / "start.extxyz", case / "start.extxyz")
    shutil.copy(SOURCE_EXAMPLE / "mace_kwargs.json", case / "mace_kwargs.json")
    shutil.copy(BEC_SOURCE, case / BEC_SOURCE.name)
    _field_file(case / "field.py")

    address = f"eda-mace-{implementation}-{ensemble}-{os.getpid()}"
    (case / "input.xml").write_text(_input_xml(implementation, ensemble, address))
    env = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
    server = subprocess.Popen(
        [sys.executable, str(REPO_ROOT / "bin/i-pi"), "input.xml"],
        cwd=case,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    client = None
    try:
        _wait_for_socket(address, server)
        client = subprocess.Popen(
            [
                sys.executable,
                str(DRIVER),
                "--address",
                address,
                "--template",
                str(case / "start.extxyz"),
                "--model",
                str(model),
                "--mace-kwargs",
                str(case / "mace_kwargs.json"),
                "--bec",
                str(case / BEC_SOURCE.name),
            ],
            cwd=case,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        client_output, _ = client.communicate(timeout=180)
        server_output, _ = server.communicate(timeout=30)
        if client.returncode or server.returncode:
            raise RuntimeError(
                f"{implementation}-{ensemble} failed:\n{client_output}\n{server_output}"
            )

        from ase.io import read

        positions = np.asarray(
            [a.positions for a in read(case / "simulation.positions_0.xyz", index=":")]
        )
        velocities = np.asarray(
            [
                a.get_velocities()
                for a in read(case / "simulation.velocities_0.xyz", index=":")
            ]
        )
        return positions, velocities
    finally:
        if client is not None and client.poll() is None:
            client.kill()
        if server.poll() is None:
            server.kill()
        socket_path = Path("/tmp") / f"ipi_{address}"
        socket_path.unlink(missing_ok=True)


@pytest.mark.parametrize("ensemble", ["nve", "nvt"])
def test_eda_and_ffdielectric_match_mace_polar_dynamics(
    tmp_path, mace_polar_model, ensemble
):
    """The legacy EDA and new wrapper produce identical trajectories."""

    eda = _run_case(tmp_path, "eda", ensemble, mace_polar_model)
    dielectric = _run_case(tmp_path, "ffdielectric", ensemble, mace_polar_model)
    for legacy, new in zip(eda, dielectric):
        np.testing.assert_allclose(new, legacy, rtol=1e-11, atol=1e-12)
