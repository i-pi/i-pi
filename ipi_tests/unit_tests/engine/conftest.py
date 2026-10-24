"""Shared simulation fixtures for engine tests."""

import io
import xml.etree.ElementTree as ET
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ipi.engine.simulation import Simulation
from ipi.scripting import InteractiveSimulation
from ipi.utils.softexit import softexit
from ipi_tests.unit_tests.engine.inputs import PrescribedPolynomial, harmonic_input


@pytest.fixture
def simulation(tmp_path, monkeypatch):
    """Isolate output and close force services even when an assertion fails."""
    monkeypatch.chdir(tmp_path)
    # The soft-exit monitor polls at this interval; joining it otherwise takes ~1 s.
    monkeypatch.setattr("ipi.utils.softexit.SOFTEXITLATENCY", 0.01)

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


@pytest.fixture
def reference_system():
    """Bound read-only system with a quartic potential, nonuniform path and triclinic cell."""

    def create(count=6, alpha=0.3, style="rpmd", frequencies=()):
        root = ET.fromstring(harmonic_input(nbeads=count))
        temp = root.find("system/ensemble/temperature")
        temp.attrib.clear()
        temp.text = "0.4"
        normal = ET.SubElement(root.find("system"), "normal_modes", transform="matrix")
        ET.SubElement(normal, "frequencies", style=style).text = str(list(frequencies))
        s = Simulation.load_from_xml(
            ET.tostring(root, encoding="unicode"), read_only=True
        ).syslist[0]
        a = np.array([0.35, 0.55, 0.8, 0.45, 0.65, 0.95])
        b = np.array([0.03, 0.05, 0.02, 0.04, 0.01, 0.06])
        s.forces.ff["harmonic"].driver = PrescribedPolynomial(a, b)
        q = (
            0.35 * np.sin(np.arange(count * 6).reshape(count, 6) * 0.71)
            + np.arange(6) * 0.07
        )
        p = (
            0.23 * np.cos(np.arange(count * 6).reshape(count, 6) * 0.47)
            + np.arange(6) * 0.03
        )
        masses = np.array([2.0, 3.0])
        h = np.array([[6.0, 1.0, 0.5], [0.0, 5.0, 0.7], [0.0, 0.0, 4.0]])
        s.beads._q[:] = q
        s.beads._p[:] = p
        s.beads._m[:] = masses
        s.cell._h[:] = h
        s.forces.alpha = alpha
        return s, SimpleNamespace(
            q=q,
            p=p,
            masses=masses,
            a=a,
            b=b,
            h=h,
            temperature=0.4,
            count=count,
            alpha=alpha,
        )

    return create
