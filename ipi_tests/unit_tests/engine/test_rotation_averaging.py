"""Rotational quadrature must turn a quadratic anisotropic model into its isotropic mean."""

import numpy as np
import pytest
from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import FFDebye, FFRotations
from ipi.utils.prng import Random


@pytest.mark.parametrize("grid", ["lebedev", "legendre"])
@pytest.mark.parametrize(
    "random,inversion", [(False, False), (True, False), (False, True)]
)
def test_rotation_average_matches_trace_of_hessian(grid, random, inversion):
    hessian = np.diag([1.0, 2.0, 6.0])
    member = FFDebye(name="anisotropic", H=hessian, xref=np.zeros(3), threaded=True)
    ff = FFRotations(
        ffdirect=member,
        grid_mode=grid,
        grid_order=3,
        random=random,
        inversion=inversion,
        prng=Random(seed=17),
        dopbc=False,
    )
    atoms = Atoms(1)
    atoms.q[:] = [0.3, -0.2, 0.1]
    request = ff.queue(atoms, Cell(np.eye(3) * 10))
    ff.poll()
    assert request["status"] == "Done"
    energy, force, virial, extras = request["result"]
    k = np.trace(hessian) / 3
    assert energy == pytest.approx(k * np.dot(atoms.q, atoms.q) / 2, abs=1e-12)
    np.testing.assert_allclose(force, -k * atoms.q, atol=1e-12)
    np.testing.assert_array_equal(virial, np.zeros((3, 3)))
    assert not member.requests
    assert not ff.requests
