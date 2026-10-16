"""Cavity forces must differentiate the total light-matter energy."""

import numpy as np
import pytest
from ipi.engine.forcefields import PhotonDriver


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
