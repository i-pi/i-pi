"""Fixtures shared by the unit tests."""

import io
from contextlib import contextmanager

import pytest

import ipi.utils.array_backend as ab
from ipi.scripting import InteractiveSimulation
from ipi.utils.softexit import softexit
from ipi_tests import test_tools


@pytest.fixture
def select_backend():
    """Selects an array backend, skipping the test if its library is missing,
    and restores the numpy backend after the test."""

    def select(name, device="cpu", jit=False):
        if name != "numpy":
            pytest.importorskip(name)
        ab.set_backend(name, device=device, jit=jit)

    yield select
    ab.set_backend("numpy")


@pytest.fixture
def numpy_draws(monkeypatch):
    """Makes the torch and JAX backends draw the same random numbers as numpy."""

    from ipi.utils.prng import Random

    for name, method in test_tools.numpy_draws().items():
        monkeypatch.setattr(Random, name, method)


@pytest.fixture
def simulation(tmp_path, monkeypatch, select_backend):
    """Creates simulations in a temporary folder, as context managers that stop
    them (and their force services) at the end of the block, also when an
    assertion fails. `simulation(xml, **options)` takes the options of
    InteractiveSimulation (backend, device, jit): without them, simulations use
    the backend of --ipi-backend, unless the test selects another. Several
    simulations can be active at once."""

    select_backend(test_tools.ipi_backend, test_tools.ipi_device)
    monkeypatch.chdir(tmp_path)
    softexit.reset()
    active = []

    @contextmanager
    def create(xml, **options):
        sim = None
        try:
            sim = InteractiveSimulation(io.StringIO(xml), **options)
            active.append(sim)
            yield sim
        finally:
            if sim is not None:
                sim.stop()
                active.remove(sim)
            if not active:
                softexit.reset()

    try:
        yield create
    finally:
        softexit.reset()
