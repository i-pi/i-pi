"""Rotational quadrature must turn a quadratic anisotropic model into its isotropic mean."""

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import FFDebye, FFRotations, ForceField, ForceRequest
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
    atoms._q[:] = [0.3, -0.2, 0.1]
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

    member = ForceField(name="zero", threaded=True)
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


@pytest.mark.parametrize(
    "extras,expected",
    [
        ([{"label": "first"}, {"label": "second"}], {"label": ["first", "second"]}),
        (["first", "second"], {"raw": "[ first, second ]"}),
    ],
)
def test_rotation_average_preserves_extras_from_each_orientation(extras, expected):

    member = ForceField()
    ff = FFRotations(ffdirect=member, grid_order=1, prng=Random(seed=19))
    requests = [
        ForceRequest(result=[0.0, np.zeros(3), np.zeros((3, 3)), extra])
        for extra in extras
    ]
    request = ForceRequest(
        pos=np.zeros(3), ff_handles=requests, rots=[(np.eye(3), 0.5), (-np.eye(3), 0.5)]
    )
    ff.gather(request)
    for key, value in expected.items():
        assert request["result"][3][key] == value
