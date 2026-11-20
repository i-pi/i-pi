"""Rotational quadrature must turn a quadratic anisotropic model into its isotropic mean."""

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import FFDebye, FFRotations, ForceField, ForceRequest
from ipi.utils.mathtools import random_rotation
from ipi.utils.prng import Random


@pytest.mark.parametrize("grid", ["lebedev", "legendre"])
@pytest.mark.parametrize(
    "random,inversion", [(False, False), (True, False), (False, True)]
)
def test_rotation_average_matches_trace_of_hessian(grid, random, inversion):
    hessian = np.diag([1.0, 2.0, 6.0])
    member = FFDebye(name="anisotropic", H=hessian, xref=np.zeros(3))
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
    atoms.q_[:] = [0.3, -0.2, 0.1]
    request = ff.queue(atoms, Cell(np.eye(3) * 10))
    member.poll()
    ff.poll()
    assert request["status"] == "Done"
    energy, force, virial, extras = request["result"]
    k = np.trace(hessian) / 3
    assert energy == pytest.approx(k * np.dot(atoms.q, atoms.q) / 2, abs=1e-12)
    np.testing.assert_allclose(force, -k * atoms.q, atol=1e-12)
    np.testing.assert_array_equal(virial, np.zeros((3, 3)))
    assert not member.requests
    assert not ff.requests


def test_rotation_average_completes_all_ready_requests_in_one_poll():
    """A single poll completes all the requests whose member results are ready"""
    member = FFDebye(name="debye", H=np.eye(3), xref=np.zeros(3))
    ff = FFRotations(
        ffdirect=member,
        grid_mode="lebedev",
        grid_order=3,
        prng=Random(seed=17),
        dopbc=False,
    )
    requests = []
    for i in range(4):
        atoms = Atoms(1)
        atoms.q_[:] = [0.1 * (i + 1), 0.0, 0.0]
        requests.append(ff.queue(atoms, Cell(np.eye(3) * 10)))
    member.poll()
    ff.poll()
    assert [r["status"] for r in requests] == ["Done"] * 4
    assert not ff.requests


@pytest.mark.parametrize(
    "options,message",
    [
        ({}, "Must specify"),
        ({"ffsocket": True, "ffdirect": True}, "Cannot specify both"),
        ({"ffdirect": True, "grid_mode": "invalid"}, "Invalid quadrature"),
    ],
)
def test_rotation_average_rejects_ambiguous_models(options, message):

    options = {
        key: ForceField(name=key) if value is True else value
        for key, value in options.items()
    }
    with pytest.raises(ValueError, match=message):
        FFRotations(**options)


def test_rotation_average_waits_for_all_member_results():
    """Partial results must not be reported as a completed rotation average."""

    member = ForceField(name="zero")
    ff = FFRotations(ffsocket=member, grid_order=1, prng=Random(seed=19))
    ff.bind(None)
    request = ff.queue(Atoms(1), Cell(np.eye(3)))
    ff.poll()
    assert request["status"] == "Running"
    assert request["result"] is None
    member.poll()
    ff.poll()
    assert request["status"] == "Done"
    np.testing.assert_array_equal(request["result"][1], np.zeros(3))
    assert not member.requests


def test_rotation_average_polled_while_queued(monkeypatch):
    """A poll may see a request while it is being queued, possibly with all
    its member results: it must complete it only once queued."""

    queue = ForceField.queue

    def queue_and_poll(self, *args, **kwargs):
        request = queue(self, *args, **kwargs)
        if isinstance(self, FFRotations):
            self.poll()
        return request

    monkeypatch.setattr(ForceField, "queue", queue_and_poll)
    member = FFDebye(name="isotropic", H=np.eye(3), xref=np.zeros(3))
    ff = FFRotations(ffdirect=member, random=True, prng=Random(seed=23), dopbc=False)
    request = ff.queue(Atoms(1), Cell(np.eye(3) * 10))
    member.poll()
    ff.poll()
    assert request["status"] == "Done"
    assert not ff.requests


@pytest.mark.parametrize(
    "extras,expected",
    [
        ([{"label": "first"}, {"label": "second"}], {"label": ["first", "second"]}),
        (["first", "second"], {"raw": "[ first, second ]"}),
    ],
)
def test_rotation_average_preserves_extras_from_each_orientation(extras, expected):

    member = ForceField()
    ff = FFRotations(
        ffdirect=member, grid_order=1, inversion=True, prng=Random(seed=19)
    )
    requests = [
        ForceRequest(result=[0.0, np.zeros(3), np.zeros((3, 3)), extra])
        for extra in extras
    ]
    request = ForceRequest(
        pos=np.zeros(3), ff_handles=requests, rots=np.stack([np.eye(3), -np.eye(3)])
    )
    ff.gather(request)
    for key, value in expected.items():
        assert request["result"][3][key] == value


def test_random_rotations_are_uniform():
    """Random rotations are orthogonal, improper half of the time, and their
    proper part is uniform: its elements average to 0, their squares to 1/3"""

    prng = Random(seed=11)
    R = np.array([random_rotation(prng) for _ in range(4000)])
    identity = np.broadcast_to(np.eye(3), R.shape)
    np.testing.assert_allclose(R @ R.transpose(0, 2, 1), identity, atol=1e-12)
    det = np.linalg.det(R)
    assert abs((det < 0).mean() - 0.5) < 0.03
    proper = R * det[:, None, None]
    assert np.abs(proper.mean(axis=0)).max() < 0.04
    np.testing.assert_allclose(
        (proper**2).mean(axis=0), np.full((3, 3), 1 / 3), atol=0.03
    )


@pytest.mark.parametrize("inversion", [False, True])
def test_rotation_average_of_backend_arrays(select_backend, inversion):
    """With the torch backend, the averages are torch arrays, also when the
    rotated structures are evaluated with numpy"""

    select_backend("torch")
    import torch

    from ipi.utils.array_backend import xp

    hessian = np.diag([1.0, 2.0, 6.0])
    member = FFDebye(name="anisotropic", H=hessian, xref=np.zeros(3))
    ff = FFRotations(
        ffdirect=member,
        grid_order=3,
        random=True,
        inversion=inversion,
        prng=Random(seed=17),
        dopbc=False,
    )
    atoms = Atoms(1)
    atoms.q_[:] = xp.asarray([0.3, -0.2, 0.1])
    request = ff.queue(atoms, Cell(xp.eye(3) * 10))
    member.poll()
    ff.poll()
    energy, force, virial, extras = request["result"]
    assert all(torch.is_tensor(x) for x in (energy, force, virial))
    k = np.trace(hessian) / 3
    q = atoms.q.numpy()
    assert float(energy) == pytest.approx(k * np.dot(q, q) / 2, abs=1e-12)
    np.testing.assert_allclose(force.numpy(), -k * q, atol=1e-12)
    assert extras["o3grid_pots"].shape == (18 * (2 if inversion else 1),)
