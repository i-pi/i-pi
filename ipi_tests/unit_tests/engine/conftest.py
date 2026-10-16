"""Shared isolated simulation fixture for numerical reference tests."""

import io
from contextlib import contextmanager
import pytest
from ipi.scripting import InteractiveSimulation
from ipi.utils.softexit import softexit


@pytest.fixture
def simulation(tmp_path, monkeypatch):
    """Isolate output and close force services even when an assertion fails."""
    monkeypatch.chdir(tmp_path)

    @contextmanager
    def create(xml):
        softexit.reset()
        sim = None
        try:
            sim = InteractiveSimulation(io.StringIO(xml))
            yield sim
        finally:
            if sim is not None:
                sim.stop()
            softexit.reset()

    try:
        yield create
    finally:
        softexit.reset()
