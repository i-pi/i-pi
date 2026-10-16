"""Built-in force providers checked against independent energies and derivatives."""

import numpy as np
import pytest
from ipi.engine.forcefields import FFLennardJones, FFDebye, FFdmd, ForceRequest


def evaluate(ff, q, h):
    request = ForceRequest(pos=q.ravel().copy(), cell=(h, np.linalg.inv(h)))
    ff.evaluate(request)
    assert request["status"] == "Done" and request._event_done.is_set()
    return request["result"]


@pytest.mark.parametrize("kind", ["lj", "debye", "dmd"])
def test_forcefield_gradient_and_translation_properties(kind):
    q = np.array([[0.1, 0.2, 0.3], [1.4, 0.3, 0.2], [0.3, 1.7, -0.1]])
    h = np.eye(3) * 10
    if kind == "lj":
        ff = FFLennardJones(
            pars={"eps": 0.7, "sigma": 1.1}, dopbc=False, threaded=False
        )
    elif kind == "debye":
        ff = FFDebye(H=np.diag(np.linspace(0.5, 1.3, 9)), xref=np.zeros(9), vref=0.2)
    else:
        ff = FFdmd(coupling=np.array([0.2, 0.3, 0.4]), freq=0.5, dtdmd=0.2, dmdstep=3)
    energy, force, vir, *_ = evaluate(ff, q, h)
    for i in range(q.size):
        delta = np.zeros_like(q)
        delta.flat[i] = 1e-6
        expected = -(evaluate(ff, q + delta, h)[0] - evaluate(ff, q - delta, h)[0]) / (
            2e-6
        )
        assert force[i] == pytest.approx(expected, abs=1e-8)
    if kind != "debye":
        np.testing.assert_allclose(force.reshape(-1, 3).sum(axis=0), 0, atol=1e-12)
        assert evaluate(ff, q + 2, h)[0] == pytest.approx(energy)
    if kind == "dmd":
        eps = 1e-6
        dv = (
            evaluate(ff, q * (1 + eps), h * (1 + eps))[0]
            - evaluate(ff, q * (1 - eps), h * (1 - eps))[0]
        ) / (2 * eps)
        assert np.trace(vir) == pytest.approx(-dv, abs=1e-9)
        ff.dmd_update()
        assert evaluate(ff, q, h)[0] / energy == pytest.approx(
            np.sin(0.4) / np.sin(0.3)
        )


@pytest.mark.parametrize("kwargs", [{}, {"H": np.eye(3)}])
def test_debye_requires_hessian_and_reference(kwargs):
    with pytest.raises(ValueError):
        FFDebye(**kwargs)


def test_lennard_jones_minimum_and_periodic_rejection():
    with pytest.raises(ValueError, match="Periodic"):
        FFLennardJones(pars={"eps": 1, "sigma": 1})
    ff = FFLennardJones(pars={"eps": 0.7, "sigma": 1.1}, dopbc=False)
    q = np.array([[0.0, 0, 0], [2 ** (1 / 6) * 1.1, 0, 0]])
    e, f, *_ = evaluate(ff, q, np.eye(3) * 10)
    assert e == pytest.approx(-0.7)
    np.testing.assert_allclose(f, 0, atol=1e-12)


@pytest.mark.parametrize("baseline", [False, True])
def test_committee_force_differentiates_combined_energy(baseline):
    from ipi.engine.forcefields import FFCommittee
    from ipi.engine.atoms import Atoms
    from ipi.engine.cell import Cell

    members = [
        FFDebye(name=str(k), H=np.eye(3) * k, xref=np.zeros(3)) for k in [1.0, 2.0, 3.0]
    ]
    ff = FFCommittee(
        fflist=members,
        baseline_name="1.0" if baseline else "",
        baseline_uncertainty=0.4 if baseline else -1.0,
        alpha=1.3,
    )
    cell = Cell(np.eye(3) * 10)
    atoms = Atoms(1)
    atoms.m[:] = 1

    def compute(q):
        atoms.q[:] = q
        request = ff.queue(atoms, cell)
        for member in members:
            member.poll()
        ff.poll()
        assert request["status"] == "Done"
        result = request["result"]
        ff.release(request)
        return result

    q = np.array([0.2, -0.3, 0.4])
    energy, force, _, extras = compute(q)
    for i in range(3):
        delta = np.zeros(3)
        delta[i] = 1e-6
        assert force[i] == pytest.approx(
            -(compute(q + delta)[0] - compute(q - delta)[0]) / (2e-6), abs=1e-9
        )
    if not baseline:
        assert energy == pytest.approx(np.dot(q, q))
        np.testing.assert_allclose(force, -2 * q)
    assert extras["committee_uncertainty"] > 0
    assert all(len(member.requests) == 0 for member in members)
