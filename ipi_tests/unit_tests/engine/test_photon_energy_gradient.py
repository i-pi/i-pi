"""Cavity forces must differentiate the total light-matter energy."""

from copy import deepcopy

import numpy as np
import pytest

from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.engine.forcefields import FFCavPhSocket, ForceField, PhotonDriver
from ipi.utils.softexit import Softexit


@pytest.mark.parametrize("representation", ["loose", "dense"])
@pytest.mark.parametrize("born_charges", [False, True])
def test_photon_and_nuclear_forces_are_energy_gradients(representation, born_charges):
    driver = PhotonDriver(E0=0.3, omega_c=0.7, ph_rep=representation)
    nuclei = np.array([[0.2, -0.3, 0.4], [-0.1, 0.5, 0.2]])
    charges = np.array([1.0, -1.0])
    photons = np.linspace(-0.2, 0.4, driver.n_photon_3)

    def energy(q, p):
        driver.split_atom_ph_coord(np.concatenate([q.ravel(), p]))
        dipole = charges @ q
        return driver.get_ph_energy(dipole[:1], dipole[1:2])

    energy(nuclei, photons)
    dipole = charges @ nuclei
    fp = driver.get_ph_forces(dipole[:1], dipole[1:2])
    kwargs = (
        {"dipole_der": np.array([c * np.eye(3) for c in charges]).ravel()}
        if born_charges
        else {"charge_array_bath": charges}
    )
    fn = np.array(driver.get_nuc_cav_forces(dipole[:1], dipole[1:2], **kwargs)).T
    eps = 1e-6
    for i in range(photons.size):
        dq = np.zeros_like(photons)
        dq[i] = eps
        assert fp[i] == pytest.approx(
            -(energy(nuclei, photons + dq) - energy(nuclei, photons - dq)) / (2 * eps),
            abs=1e-10,
        )
    for i in range(nuclei.size):
        dq = np.zeros_like(nuclei)
        dq.flat[i] = eps
        assert fn.flat[i] == pytest.approx(
            -(energy(nuclei + dq, photons) - energy(nuclei - dq, photons)) / (2 * eps),
            abs=1e-10,
        )


def test_disabled_photons_leave_all_coordinates_atomic():
    q = np.arange(12.0)
    atoms, photons = PhotonDriver(apply_photon=False).split_atom_ph_coord(q)
    np.testing.assert_array_equal(atoms, q)
    assert photons.size == 0


@pytest.mark.parametrize(
    "extras,expected",
    [
        ([], {"raw": ""}),
        ([{"dipole": [1, 2, 3]}], {"dipole": [1, 2, 3], "raw": ""}),
        (["legacy"], {"raw": "legacy"}),
        (["legacy", 2], {"raw": ["legacy", "2"]}),
        (
            [{"raw": "a", "energy": 1}, {"energy": 2}],
            {"raw": ["a", ""], "energy": [1, 2]},
        ),
        ([{"dipole": [1, 2, 3]}, {}], {"dipole": [[1, 2, 3], None], "raw": ["", ""]}),
        ([{"data": [1]}, {"data": [2, 3]}], {"data": [[1], [2, 3]], "raw": ["", ""]}),
        ([{"label": "a"}, {"label": "b"}], {"label": ["a", "b"], "raw": ["", ""]}),
    ],
)
def test_cavity_bath_extras_preserve_bath_order(extras, expected):

    original = deepcopy(extras)
    ff = FFCavPhSocket(apply_photon=False)
    actual = ff.combine_bath_extras(extras)
    assert actual.keys() == expected.keys()
    for key, value in expected.items():
        if isinstance(actual[key], np.ndarray):
            np.testing.assert_array_equal(actual[key], value)
        else:
            assert actual[key] == value
    assert extras == original


def test_cavity_dipoles_keep_independent_baths_separate():

    ff = FFCavPhSocket(apply_photon=False)
    positions = np.array([[1, 2, 3], [0, 1, 2], [4, 5, 6], [2, 1, 0]])
    dipoles = ff.calc_dipole_xyz_mm(positions.ravel(), 2, np.array([1.0, -1.0]))
    np.testing.assert_array_equal(np.array(dipoles).T, [[1, 1, 1], [2, 4, 6]])


@pytest.mark.parametrize("representation", ["loose", "dense"])
@pytest.mark.parametrize(
    "apply,evaluate,surface",
    [
        (False, False, False),
        (True, False, False),
        (True, True, False),
        (True, True, True),
    ],
)
def test_cavity_request_force_matches_energy(representation, apply, evaluate, surface):
    """Use the real zero-potential evaluator to isolate cavity contributions.

    ForceField.poll supplies molecular results in process; socket transport is
    tested separately. Queue splitting and photon assembly remain unmodified.
    """

    class PointChargeDriver(ForceField):
        """Analytic dipole surface for two fixed charges, with zero nuclear PES."""

        def poll(self):
            super().poll()
            for request in self.requests:
                q = request["pos"].reshape(-1, 3)
                request["result"][3].update(
                    dipole=q[0] - q[1],
                    dipole_derivative=np.array([np.eye(3), -np.eye(3)]).ravel(),
                )

    transport = PointChargeDriver() if surface else ForceField()
    ff = FFCavPhSocket(
        interface=transport,
        threaded=False,
        apply_photon=apply,
        evaluate_photon=evaluate,
        dipole_surface=surface,
        ph_rep=representation,
        charge_array=np.array([1.0, -1.0]),
        E0=0.3,
        omega_c=0.7,
        pars={"model": "zero"},
        dopbc=True,
    )
    count = 2 + (ff.ph.n_photon if apply else 0)
    atoms = Atoms(count)
    q = np.linspace(-0.3, 0.4, count * 3)
    cell = Cell(np.eye(3) * 20)

    def calculate(pos):
        atoms.q[:] = pos
        request = ff.queue(atoms, cell, reqid=4)
        assert request["status"] == "Done"
        assert not ff.requests
        return request["result"]

    energy, force, virial, _ = calculate(q)
    np.testing.assert_array_equal(virial, np.zeros((3, 3)))
    if not (apply and evaluate):
        assert energy == 0
        np.testing.assert_array_equal(force, np.zeros_like(q))
    else:
        # The cavity Hamiltonian is a sum of displaced harmonic oscillators.
        displaced = q[6:].copy()
        polarized = [0, 4 if representation == "loose" else 1]
        displaced[polarized] += (0.3 / 0.7**2) * (q[:2] - q[3:5])
        assert energy == pytest.approx(0.5 * 0.7**2 * np.dot(displaced, displaced))
        for i in range(q.size):
            delta = np.zeros_like(q)
            delta[i] = 1e-6
            assert force[i] == pytest.approx(
                -(calculate(q + delta)[0] - calculate(q - delta)[0]) / 2e-6, abs=1e-10
            )


@pytest.mark.parametrize("active", [[2], [0, 1, 2]])
def test_cavity_rejects_active_atoms_outside_system(active):

    ff = FFCavPhSocket(apply_photon=False, active=np.array(active))
    with pytest.raises(ValueError, match="more active atoms"):
        ff.queue(Atoms(2), Cell(np.eye(3)))
    assert not ff.requests


@pytest.mark.xfail(
    strict=True,
    raises=KeyError,
    reason="Cavity dipole validation indexes missing extras before checking their presence",
)
def test_cavity_missing_dipole_data_reports_a_validation_error():

    ff = FFCavPhSocket(interface=ForceField(), threaded=False, dipole_surface=True)
    with pytest.raises(ValueError, match="[Dd]ipole"):
        ff.queue(Atoms(1 + ff.ph.n_photon), Cell(np.eye(3) * 10))


@pytest.mark.parametrize("where", ["dipole", "queue", "cancelled-request"])
@pytest.mark.xfail(
    strict=True,
    raises=ValueError,
    reason="Cavity error handling passes the error message as an invalid soft-exit status",
)
def test_cavity_invalid_or_cancelled_request_triggers_soft_exit(where, monkeypatch):

    exit_state = Softexit()
    monkeypatch.setattr("ipi.engine.forcefields.softexit", exit_state)
    transport = ForceField()
    if where == "cancelled-request":
        monkeypatch.setattr(transport, "poll", transport.stop)
    ff = FFCavPhSocket(interface=transport, threaded=False, charge_array=np.ones(2))
    with pytest.raises(SystemExit):
        if where == "dipole":
            ff.calc_dipole_xyz_mm(np.zeros(3), 1, np.ones(2))
        else:
            ff.queue(Atoms(1 + ff.ph.n_photon), Cell(np.eye(3) * 10))
    assert exit_state.triggered
