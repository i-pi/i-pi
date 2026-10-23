"""Unit tests for the FFDirect forcefield."""

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import FFDirect

K = 2.0


def _make_atoms(q):
    atoms = Atoms(len(q) // 3)
    atoms.q[:] = np.asarray(q, dtype=float)
    return atoms


def _assert_harmonic_result(req, q):
    energy, force, virial, extra = req["result"]
    assert np.isclose(energy, 0.5 * K * (q**2).sum())
    assert np.allclose(force, -K * q)
    assert np.allclose(virial, 0.0)
    assert extra == {"raw": "nada"}


def test_ffdirect_harmonic_single_request():
    """Checks a single in-process harmonic force evaluation."""

    q = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ff = FFDirect(
        name="direct",
        pes="harmonic",
        pars={"k1": K},
        dopbc=False,
        threaded=False,
    )

    req = ff.queue(_make_atoms(q), Cell(np.eye(3) * 10.0), reqid=0)

    assert req["status"] == "Done"
    _assert_harmonic_result(req, q)


def test_ffdirect_harmonic_batched_requests():
    """Checks that a full FFDirect batch returns correct harmonic forces."""

    cell = Cell(np.eye(3) * 10.0)
    positions = [
        np.array([0.1, 0.0, 0.0]),
        np.array([0.0, 0.2, 0.0]),
        np.array([0.0, 0.0, 0.3]),
        np.array([0.4, -0.5, 0.6]),
    ]
    ff = FFDirect(
        name="direct_batch",
        pes="harmonic",
        pars={"k1": K},
        dopbc=False,
        threaded=False,
        batch_size=len(positions),
    )

    reqs = [ff.queue(_make_atoms(q), cell, reqid=i) for i, q in enumerate(positions)]

    for req, q in zip(reqs, positions):
        assert req["status"] == "Done"
        _assert_harmonic_result(req, q)


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="BUG: gas FFDirect batches use the first structure's atom count for every force array",
)
def test_ffdirect_gas_batch_preserves_different_atom_counts():
    """A shared gas provider must return one force vector per input coordinate."""
    positions = [
        np.array([0.2, -0.3, 0.1, -0.4, 0.2, 0.3]),
        np.array([0.1, 0.2, 0.3, -0.2, 0.4, -0.1, 0.5, -0.3, 0.2]),
    ]
    ff = FFDirect(pes="gas", batch_size=2, threaded=False)
    cell = Cell(np.eye(3) * 20)
    requests = [
        ff.queue(_make_atoms(q), cell, reqid=i) for i, q in enumerate(positions)
    ]

    assert len(requests) == 2
    for request, q in zip(requests, positions):
        assert request["status"] == "Done"
        assert request._event_done.is_set()
        energy, force, virial, _ = request["result"]
        assert force.shape == q.shape
        np.testing.assert_array_equal(force, np.zeros(q.shape))
        assert energy == 0.0
        np.testing.assert_array_equal(virial, np.zeros((3, 3)))
        ff.release(request)


@pytest.mark.xfail(
    strict=True,
    raises=ValueError,
    reason="BUG: harmonic FFDirect batches call np.stack on structures with different atom counts",
)
def test_ffdirect_harmonic_batch_preserves_different_atom_counts():
    """A shared harmonic provider must preserve each structure's size and results."""
    positions = [
        np.array([0.2, -0.3, 0.1, -0.4, 0.2, 0.3]),
        np.array([0.1, 0.2, 0.3, -0.2, 0.4, -0.1, 0.5, -0.3, 0.2]),
    ]
    spring = np.array([2.0, 3.0, 5.0])
    ff = FFDirect(
        pes="harmonic",
        pars={"k1": 2.0, "k2": 3.0, "k3": 5.0},
        batch_size=2,
        threaded=False,
    )
    cell = Cell(np.eye(3) * 20)
    requests = [
        ff.queue(_make_atoms(q), cell, reqid=i) for i, q in enumerate(positions)
    ]

    assert len(requests) == 2
    for request, q in zip(requests, positions):
        assert request["status"] == "Done"
        assert request._event_done.is_set()
        energy, force, virial, _ = request["result"]
        assert force.shape == q.shape
        coordinates = q.reshape(-1, 3)
        np.testing.assert_allclose(force, (-spring * coordinates).reshape(-1))
        assert energy == pytest.approx(0.5 * np.sum(spring * coordinates**2))
        np.testing.assert_array_equal(virial, np.zeros((3, 3)))
        ff.release(request)


def test_ffdirect_flushes_incomplete_batch_after_idle_polls():
    """A final short batch must finish even if no more structures arrive."""
    ff = FFDirect(pes="harmonic", pars={"k1": K}, batch_size=3, threaded=True)
    q = np.array([0.1, -0.2, 0.3])
    request = ff.queue(_make_atoms(q), Cell(np.eye(3) * 10))
    ff.poll()
    assert request["status"] == "Running"
    assert not request._event_done.is_set()
    for _ in range(20):
        if request._event_done.is_set():
            break
        ff.poll()
    assert request._event_done.is_set()
    _assert_harmonic_result(request, q)
    ff.release(request)


def test_ffdirect_requires_custom_driver_path():
    with pytest.raises(ValueError, match="pes_path"):
        FFDirect(pes="custom")


def test_ffdirect_preserves_driver_import_error(tmp_path):
    path = tmp_path / "unavailable_driver.py"
    path.write_text("raise ImportError('required model library is unavailable')\n")
    with pytest.raises(ImportError, match="required model library"):
        FFDirect(pes="custom", pes_path=str(path))
