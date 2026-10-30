"""Built-in force providers checked against energy derivatives and known limits."""

import json

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import (
    FFCommittee,
    FFDebye,
    FFdmd,
    FFLennardJones,
    ForceField,
    ForceRequest,
)
from ipi.engine.outputs import OutputMaker


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


@pytest.mark.parametrize(
    "options,message",
    [
        ({"fflist": []}, "empty list"),
        ({"ffweights": [1]}, "weights"),
        ({"baseline_name": "base"}, "simultaneously"),
        ({"baseline_uncertainty": 1}, "simultaneously"),
        (
            {"baseline_name": "base", "baseline_uncertainty": 1, "ffweights": [1, 1]},
            "weights",
        ),
    ],
)
def test_committee_rejects_inconsistent_configuration(options, message):

    defaults = dict(fflist=[ForceField(name="base"), ForceField(name="member")])
    with pytest.raises(ValueError, match=message):
        FFCommittee(**(defaults | options))


@pytest.mark.parametrize(
    "H,xref,message",
    [
        (np.eye(6), np.zeros(6), "Hessian"),
        (np.eye(3), np.zeros(6), "Reference"),
    ],
)
def test_debye_rejects_request_with_incompatible_dimensions(H, xref, message):
    ff = FFDebye(H=H, xref=xref)
    with pytest.raises(ValueError, match=message):
        evaluate(ff, np.zeros(3), np.eye(3))


@pytest.mark.parametrize("components", ["both", "force", "virial", "neither"])
def test_committee_json_ensemble_preserves_means_and_uncertainty(components):
    """A serialized two-model ensemble has known mean and sample variance."""

    extras = {"committee_pot": [1.0, 3.0], "raw": "two models"}
    if components in ("both", "force"):
        extras["committee_force"] = [[1, 2, 3], [3, 4, 5]]
    if components in ("both", "virial"):
        extras["committee_virial"] = [np.eye(3), 3 * np.eye(3)]
    member = ForceField()
    ff = FFCommittee(fflist=[member], parse_json=True, alpha=2)
    request = ForceRequest(
        pos=np.zeros(3),
        ff_handles=[
            ForceRequest(result=[2.0, np.array([2, 3, 4]), 2 * np.eye(3), extras])
        ],
    )
    ff.gather(request)
    energy, force, virial, combined = request["result"]
    assert energy == 2
    np.testing.assert_array_equal(force, [2, 3, 4])
    np.testing.assert_array_equal(virial, 2 * np.eye(3))
    np.testing.assert_array_equal(combined["committee_pot"], [0, 4])
    assert combined["committee_uncertainty"] == pytest.approx(np.sqrt(8))
    assert ("committee_force" in combined) == (components in ("both", "force"))
    assert ("committee_virial" in combined) == (components in ("both", "virial"))
    if components in ("both", "force"):
        np.testing.assert_array_equal(
            combined["committee_force"], [[0, 1, 2], [4, 5, 6]]
        )


@pytest.mark.parametrize(
    "problem,message",
    [
        ("force", "all.*components"),
        ("virial", "all.*components"),
        ("baseline", "force ensemble"),
        ("collision", "Name clash"),
    ],
)
def test_committee_rejects_ambiguous_member_results(problem, message):

    def result(extra):
        return ForceRequest(result=[0.0, np.zeros(3), np.zeros((3, 3)), extra])

    extras = {"committee_pot": [1.0, 3.0]}
    if problem == "virial":
        extras["committee_force"] = [np.zeros(3), np.zeros(3)]
    handles = [result(extras), result({"raw": ""})]
    options = {}
    if problem == "baseline":
        options = dict(baseline_name="base", baseline_uncertainty=1.0)
        handles.reverse()
    elif problem == "collision":
        handles = [result({}), result({"pot": 2})]
    ff = FFCommittee(
        fflist=[ForceField(name="base"), ForceField(name="model")],
        parse_json=True,
        **options
    )
    with pytest.raises(ValueError, match=message):
        ff.gather(ForceRequest(pos=np.zeros(3), ff_handles=handles))


def test_committee_saves_uncertain_structures(tmp_path):

    members = [FFDebye(H=np.eye(3) * k, xref=np.zeros(3)) for k in [1.0, 3.0]]
    ff = FFCommittee(fflist=members, active_thresh=0.1)
    with pytest.raises(ValueError, match="output file"):
        ff.bind(OutputMaker())
    ff.active_out = str(tmp_path / "uncertain.json")
    ff.bind(OutputMaker())
    try:
        atoms = Atoms(1)
        atoms.q[:] = [1, 0, 0]
        request = ff.queue(atoms, Cell(np.eye(3) * 10))
        ff.poll()
        ff.active_file.force_flush()
        saved = json.loads((tmp_path / "uncertain.json").read_text())
        assert saved["position"] == [1, 0, 0]
        assert saved["cell"] == (np.eye(3) * 10).ravel().tolist()
        assert saved["uncertainty"] == pytest.approx(np.sqrt(0.5))
        assert request._event_done.is_set()
        assert all(not member.requests for member in members)
    finally:
        ff.active_file.close_stream()


@pytest.mark.parametrize("missing", ["coupling", "freq", "dtdmd"])
def test_dmd_requires_couplings_frequency_and_timestep(missing):
    options = dict(coupling=np.array([1.0]), freq=0.5, dtdmd=0.1)
    options[missing] = None
    with pytest.raises(ValueError, match="Must provide"):
        FFdmd(**options)


def test_dmd_rejects_wrong_number_of_pair_couplings():
    ff = FFdmd(coupling=np.ones(2), freq=0.5, dtdmd=0.1)
    with pytest.raises(ValueError, match="Coupling matrix size"):
        evaluate(ff, np.zeros((3, 3)), np.eye(3))


def test_committee_waits_for_every_member():

    members = [ForceField(threaded=True), ForceField(threaded=True)]
    ff = FFCommittee(fflist=members)
    request = ff.queue(Atoms(1), Cell(np.eye(3)))
    members[0].poll()
    ff.poll()
    assert request["result"] is None
    assert not request._event_done.is_set()
    members[1].poll()
    ff.poll()
    assert request._event_done.is_set()
    assert request["result"][0] == 0
    assert request["result"][3]["committee_uncertainty"] == 0
    assert all(not member.requests for member in members)
