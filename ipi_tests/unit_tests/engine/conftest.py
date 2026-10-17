"""Shared isolated simulation fixture for numerical reference tests."""

import io
from contextlib import contextmanager
from pathlib import Path
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
            # Stop monitors before reset replaces the shared loop controls.
            for thread, loop_control in softexit.tlist:
                if loop_control is not None:
                    loop_control[0] = False
            for thread, _ in softexit.tlist:
                thread.join()
            softexit.reset()

    try:
        yield create
    finally:
        softexit.reset()


@pytest.fixture
def restart_states(simulation):
    """Return snapshots after continuous and interrupted runs of equal length."""

    def compare(xml, snapshot, steps=(2, 4)):
        before, after = steps
        with simulation(xml) as sim:
            sim.run(steps=before + after, write_outputs=False)
            continuous = snapshot(sim)
        with simulation(xml) as sim:
            sim.run(steps=before, write_outputs=False)
        with simulation(Path("RESTART").read_text()) as sim:
            assert sim.step == before
            sim.run(steps=after, write_outputs=False)
            restarted = snapshot(sim)
        return continuous, restarted

    return compare
