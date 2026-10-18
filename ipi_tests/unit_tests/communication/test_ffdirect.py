"""Unit tests for the FFDirect forcefield."""

import numpy as np

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


def test_ffdirect_flushes_incomplete_batch_after_idle_polls():
    """A final short batch must finish even if no more structures arrive."""
    ff = FFDirect(pes="harmonic", pars={"k1": K}, batch_size=3, threaded=True)
    q = np.array([0.1, -0.2, 0.3])
    request = ff.queue(_make_atoms(q), Cell(np.eye(3) * 10))
    ff.poll()
    assert request["status"] == "Running"
    assert not request._event_done.is_set()
    for _ in range(3):
        ff.poll()
    assert request._event_done.is_set()
    _assert_harmonic_result(request, q)
    assert not ff.request_batch
    ff.release(request)


def test_ffdirect_requires_custom_driver_path():
    import pytest

    with pytest.raises(ValueError, match="pes_path"):
        FFDirect(pes="custom")


def test_ffdirect_preserves_driver_import_error(tmp_path):
    import pytest

    path = tmp_path / "unavailable_driver.py"
    path.write_text("raise ImportError('required model library is unavailable')\n")
    with pytest.raises(ImportError, match="required model library"):
        FFDirect(pes="custom", pes_path=str(path))
