"""Checks that the numpy, torch and JAX backends give the same dynamics, and
that an MD step does not synchronize the host with a CUDA device."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import xml.etree.ElementTree as ET

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi_tests.unit_tests.engine.inputs import (
    CASES,
    HARMONIC,
    NO_FORCES,
    harmonic_input,
    write_input,
)

try:
    import torch

    cuda_available = torch.cuda.is_available()
except ImportError:
    cuda_available = False

# the numpy trajectories, computed once for each input
_REFERENCES = {}


def trajectory(simulation, xml, backend, device, nsteps):
    """Runs nsteps and returns positions, momenta, cell and conserved energy."""

    frames = []
    with simulation(xml, backend=backend, device=device) as sim:
        for step in range(nsteps):
            sim.run_step(step)
            system = sim.syslist[0]
            frames.append(
                [
                    np.array(ab.to_numpy(system.beads.q)),
                    np.array(ab.to_numpy(system.beads.p)),
                    np.array(ab.to_numpy(system.cell.h)),
                    float(sim.properties("conserved")),
                ]
            )
    return frames


def reference(simulation, xml, nsteps):
    """The numpy trajectory of an input, computed once for all the backends."""

    if (xml, nsteps) not in _REFERENCES:
        _REFERENCES[xml, nsteps] = trajectory(simulation, xml, "numpy", "cpu", nsteps)
    return _REFERENCES[xml, nsteps]


BACKENDS = [("torch", "cpu"), ("jax", "cpu")] + (
    [("torch", "cuda")] if cuda_available else []
)


@pytest.mark.parametrize("backend,device", BACKENDS)
@pytest.mark.parametrize("case", list(CASES))
def test_backend_equivalence(case, backend, device, simulation, numpy_draws, request):
    """torch and JAX follow the numpy trajectory when they draw the same
    random numbers"""

    pytest.importorskip(backend)
    if backend == "torch":
        request.getfixturevalue("strict_torch")
    if backend == "jax" and case in ("nvt_gle", "npt_gle_barostat"):
        pytest.xfail("the GLE thermostat modifies its auxiliary momenta in place")
    xml = write_input(case, HARMONIC)
    expected = reference(simulation, xml, 10)
    result = trajectory(simulation, xml, backend, device, 10)
    for ref, res in zip(expected, result):
        for a, b in zip(ref, res):
            assert np.allclose(a, b, rtol=1e-10, atol=1e-12)


# estimators computed from the positions and the forces of the beads
ESTIMATORS = [
    ("properties", "kinetic_cv"),
    ("properties", "kinetic_cv(1)"),
    ("properties", "kinetic_tens"),
    ("properties", "kinetic_tens(2)"),
    ("properties", "kinetic_ij(1;2)"),
    ("properties", "kinetic_ij(3;3)"),
    ("trajs", "kinetic_cv"),
    ("trajs", "kinetic_od"),
    ("trajs", "r_gyration"),
]


@pytest.mark.parametrize("backend,device", BACKENDS)
@pytest.mark.parametrize("case", ["nvt_langevin", "pimd_nvt"])
def test_estimators_backend_equivalence(
    case, backend, device, simulation, numpy_draws, request
):
    """torch and JAX compute the same kinetic energy estimators and gyration
    radii as numpy"""

    pytest.importorskip(backend)
    if backend == "torch":
        request.getfixturevalue("strict_torch")
    xml = write_input(case, HARMONIC)
    values = []
    for name, dev in [("numpy", "cpu"), (backend, device)]:
        with simulation(xml, backend=name, device=dev) as sim:
            sim.run_step(0)
            system = sim.syslist[0]
            values.append(
                [
                    np.array(ab.to_numpy(getattr(system, kind)[key][0]))
                    for kind, key in ESTIMATORS
                ]
            )
    for estimator, a, b in zip(ESTIMATORS, *values):
        assert np.allclose(a, b, rtol=1e-10, atol=1e-12), estimator


@pytest.mark.parametrize("multi", [False, True])
def test_jax_refuses_unsupported_motion(multi, simulation):
    """A motion that does not work with JAX, also within a multi motion, stops
    with an error when the simulation is set up"""

    pytest.importorskip("jax")
    root = ET.fromstring(harmonic_input())
    system = root.find("system")
    dynamics = system.find("motion")
    system.remove(dynamics)
    motion = ET.SubElement(system, "motion", mode="minimize")
    ET.SubElement(motion, "optimizer", mode="sd")
    if multi:
        system.remove(motion)
        ET.SubElement(system, "motion", mode="multi").extend([dynamics, motion])
    with pytest.raises(ValueError, match="GeopMotion cannot run with the JAX"):
        with simulation(ET.tostring(root, encoding="unicode"), backend="jax"):
            pass


@pytest.mark.skipif(not cuda_available, reason="CUDA not available")
@pytest.mark.parametrize(
    "case",
    [
        "nve",
        "nvt_langevin",
        "nvt_svr",
        "pimd_nvt",
        "pimd_pile_g",
        "pimd_npt",
        "npt_flexible",
    ],
)
def test_cuda_step_is_sync_free(case, simulation, device_random):
    """An MD step without force evaluation does not wait for the device"""

    xml = write_input(case, NO_FORCES)
    with simulation(xml, backend="torch", device="cuda") as sim:
        try:
            # the first steps compute and cache the integrator constants
            for step in range(2):
                sim.run_step(step)
            torch.cuda.set_sync_debug_mode("error")
            for step in range(2, 5):
                sim.run_step(step)
        finally:
            torch.cuda.set_sync_debug_mode("default")
