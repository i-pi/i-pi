"""Tests for removal of rigid motion from dynamical matrices."""

import numpy as np

from ipi.engine.beads import Beads
from ipi.utils.phonontools import apply_asr


def _nonlinear_beads():
    beads = Beads(3, 1)
    beads.m[:] = [1.0, 1.0, 1.0]
    beads.q[:] = [[-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]]
    return beads


def test_crystal_asr_projector_removes_translations():
    """Checks removal of the three mass-weighted translations."""

    beads = _nonlinear_beads()
    identity = np.eye(3 * beads.natoms)
    projector = apply_asr("crystal", identity, beads, return_trans_matrix=True)
    translations = np.zeros((3, 3 * beads.natoms))
    square_root_masses = np.sqrt(beads.m)
    for direction in range(3):
        translations[direction, direction::3] = square_root_masses

    assert np.allclose(projector, projector.T)
    assert np.allclose(projector @ projector, projector)
    assert np.allclose(projector @ translations.T, 0.0, atol=1.0e-12)
    assert np.allclose(apply_asr("crystal", identity, beads), projector)
