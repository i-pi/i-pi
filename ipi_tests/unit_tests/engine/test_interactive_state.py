"""Python/ASE access preserves physical units and invalidates dependent forces."""

import copy
import xml.etree.ElementTree as ET

import numpy as np
import pytest
import ipi.utils.array_backend as ab
from ipi.utils.units import unit_to_user
from ipi_tests.unit_tests.engine.inputs import harmonic_input


@pytest.mark.parametrize("systems,beads", [(1, 1), (1, 4), (2, 1), (2, 4)])
def test_ase_roundtrip_and_coordinate_update(simulation, systems, beads):
    root = ET.fromstring(harmonic_input(nbeads=beads))
    if systems == 2:
        root.append(copy.deepcopy(root.find("system")))
    with simulation(ET.tostring(root, encoding="unicode")) as sim:
        snapshots = [
            (s.beads.q.copy(), s.beads.p.copy(), s.cell.h.copy()) for s in sim.syslist
        ]
        structures = sim.get_structures()
        sim.set_structures(structures)
        for s, (q, p, h) in zip(sim.syslist, snapshots):
            np.testing.assert_allclose(s.beads.q, q)
            np.testing.assert_allclose(s.beads.p, p)
            np.testing.assert_allclose(s.cell.h, h)
        nested = [structures] if systems == 1 else structures
        for structure in nested:
            for atom in ([structure] if beads == 1 else structure):
                q = atom.positions / unit_to_user("length", "ase", 1)
                assert atom.get_potential_energy() == pytest.approx(
                    0.25 * np.sum(q * q) * unit_to_user("energy", "ase", 1)
                )
                np.testing.assert_allclose(
                    atom.get_forces(), -0.5 * q * unit_to_user("force", "ase", 1)
                )
                atom.positions *= 2
        sim.set_structures(structures)
        for s, (q, _, _) in zip(sim.syslist, snapshots):
            np.testing.assert_allclose(s.beads.q, 2 * q)
            np.testing.assert_allclose(s.forces.f, -q)
        values = sim.properties("potential")
        expected = [
            0.25 * np.sum(s.beads.q**2) / beads * unit_to_user("energy", "ase", 1)
            for s in sim.syslist
        ]
        np.testing.assert_allclose(values, expected[0] if systems == 1 else expected)


@pytest.mark.parametrize("temperature", [None, 0, 600])
def test_thermalization_updates_energy_accounting(simulation, temperature):
    with simulation(harmonic_input(nbeads=4)) as sim:
        system = sim.syslist[0]
        conserved = system.ensemble.econs
        sim.thermalize_momenta(temperature)
        assert system.ensemble.econs == pytest.approx(conserved)
        if temperature == 0:
            np.testing.assert_array_equal(system.beads.p, np.zeros((4, 6)))
        else:
            assert np.isfinite(system.beads.p).all()
            assert system.nm.kin > 0


def test_thermalization_with_torch_arrays(simulation, strict_torch):
    """Thermalizing works with the arrays of the torch backend, without implicit
    conversions to numpy (which fail for arrays on a GPU)"""
    ab.set_backend("torch")
    try:
        with simulation(harmonic_input(nbeads=4)) as sim:
            system = sim.syslist[0]
            conserved = float(system.ensemble.econs)
            sim.thermalize_momenta(600)
            assert float(system.ensemble.econs) == pytest.approx(conserved)
            assert float(system.nm.kin) > 0
    finally:
        ab.set_backend("numpy")


def test_custom_motion_replaces_step_and_can_be_restored(simulation):
    with simulation(harmonic_input()) as sim:
        original = sim.get_motion_step()
        q = sim.syslist[0].beads.q.copy()

        def translate(motion, step=None):
            motion.beads.q_[:] += 0.1

        sim.set_motion_step(translate)
        sim.run(steps=2, write_outputs=False)
        np.testing.assert_allclose(sim.syslist[0].beads.q, q + 0.2)
        sim.set_motion_step(original.__func__)
        assert sim.get_motion_step().__func__ is original.__func__
