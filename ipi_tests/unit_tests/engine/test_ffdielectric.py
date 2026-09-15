"""Tests for FFDielectric field application checks."""

import numpy as np
import pytest

from ipi.engine.forcefields import FFDielectric, ForceField, ForceRequest


def _client_dielectric():
    """Create the minimum FFDielectric object needed for result validation."""
    wrapped = ForceField()
    dielectric = FFDielectric(
        name="dielectric",
        where="client",
        dipole={
            "family": "electric-dipole",
            "units": "atomic_unit",
            "key": "dipole",
        },
        bec={"family": "charge", "units": "e", "key": "BEC"},
        piezo={
            "family": "electric-polarization",
            "units": "atomic_unit",
            "key": "piezoelectric",
        },
        electric_fields=[],
        electric_displacements=[],
        forcefield=wrapped,
    )
    assert not dielectric.forcefield.dopbc
    return dielectric


def test_dielectric_rejects_mixed_electric_field_and_displacement():
    """Electric-field and displacement boundary conditions are exclusive."""
    with pytest.raises(
        ValueError,
        match="cannot mix <electric_field> and <electric_displacement>",
    ):
        FFDielectric(
            name="dielectric",
            where="client",
            dipole={"family": "dipole", "units": "atomic_unit", "key": "dipole"},
            bec={"family": "born_charge", "units": "atomic_unit", "key": "BEC"},
            piezo={
                "family": "electric-polarization",
                "units": "atomic_unit",
                "key": "piezoelectric",
            },
            electric_fields=[object()],
            electric_displacements=[object()],
            forcefield=ForceField(),
        )


def test_dielectric_rejects_wrapping_in_its_wrapped_forcefield():
    """FFDielectric must reject clients configured to wrap coordinates."""
    wrapped = ForceField(dopbc=True)
    with pytest.raises(ValueError, match="does not support pbc='True'"):
        FFDielectric(
            name="dielectric",
            where="client",
            dipole={
                "family": "electric-dipole",
                "units": "atomic_unit",
                "key": "dipole",
            },
            bec={"family": "charge", "units": "e", "key": "BEC"},
            piezo={
                "family": "electric-polarization",
                "units": "atomic_unit",
                "key": "piezoelectric",
            },
            electric_fields=[],
            electric_displacements=[],
            forcefield=wrapped,
        )


def test_dielectric_expands_voigt_piezoelectric_tensor():
    """Voigt piezoelectric data is expanded with symmetric strain indices."""
    voigt = np.arange(18).reshape(3, 6)

    tensor = _client_dielectric()._coerce_piezoelectric_tensor(voigt)

    np.testing.assert_array_equal(tensor[:, 0, 0], voigt[:, 0])
    np.testing.assert_array_equal(tensor[:, 1, 1], voigt[:, 1])
    np.testing.assert_array_equal(tensor[:, 2, 2], voigt[:, 2])
    np.testing.assert_array_equal(tensor[:, 1, 2], voigt[:, 3])
    np.testing.assert_array_equal(tensor[:, 0, 2], voigt[:, 4])
    np.testing.assert_array_equal(tensor[:, 0, 1], voigt[:, 5])
    np.testing.assert_array_equal(tensor, tensor.swapaxes(1, 2))


def test_dielectric_rejects_nonsymmetric_cartesian_piezoelectric_tensor():
    """Full Cartesian piezoelectric tensors must be symmetric in strain axes."""
    tensor = np.zeros((3, 3, 3))
    tensor[0, 1, 2] = 1.0

    with pytest.raises(ValueError, match="symmetric in its two strain indices"):
        _client_dielectric()._coerce_piezoelectric_tensor(tensor)


def test_dielectric_expands_voigt_epsilon_infinity():
    """The dielectric tensor accepts six-component Voigt input."""
    epsilon = _client_dielectric()._coerce_epsilon_infinity(
        np.array([1.0, 2.0, 3.0, 0.4, 0.5, 0.6])
    )

    np.testing.assert_array_equal(
        epsilon,
        np.array([[1.0, 0.6, 0.5], [0.6, 2.0, 0.4], [0.5, 0.4, 3.0]]),
    )


def test_dielectric_rejects_nonsymmetric_epsilon_infinity():
    """A full dielectric tensor must be symmetric."""
    epsilon = np.eye(3)
    epsilon[0, 1] = 0.1

    with pytest.raises(ValueError, match="epsilon_infinity tensor must be symmetric"):
        _client_dielectric()._coerce_epsilon_infinity(epsilon)


def test_dielectric_applies_constant_displacement_ensemble():
    """Fixed D uses the screened field in the force and proper-stress terms."""
    dielectric = _client_dielectric()
    volume = 2.0
    dipole = np.array([0.1, 0.0, 0.0])
    bec = np.zeros((1, 3, 3))
    bec[0, 0, 0] = 2.0
    dfield = np.array([4.0 * np.pi * dipole[0] / volume + 2.0, 0.0, 0.0])
    request = ForceRequest(
        {
            "Dfield": dfield,
            "cell": (np.diag([volume, 1.0, 1.0]), np.eye(3)),
            "result": (
                0.0,
                np.zeros(3),
                np.zeros((3, 3)),
                {
                    "dipole": dipole,
                    "BEC": bec,
                    "piezoelectric": np.zeros((3, 3, 3)),
                    "epsilon_infinity": 2.0 * np.eye(3),
                },
            ),
        }
    )

    energy, forces, virial, _ = dielectric.fixed_D(request)

    assert energy == pytest.approx(1.0 / (2.0 * np.pi))
    np.testing.assert_array_equal(forces, [2.0, 0.0, 0.0])
    np.testing.assert_array_equal(virial, np.zeros((3, 3)))
    np.testing.assert_array_equal(dielectric.get_electric_field(), [2.0, 0.0, 0.0])


def test_dielectric_reports_displacement_for_constant_electric_field():
    """Fixed E exposes E + 4 pi mu / Omega as displacement."""
    dielectric = _client_dielectric()
    volume = 2.0
    dipole = np.array([0.1, 0.0, 0.0])
    electric_field = np.array([2.0, 0.0, 0.0])
    request = ForceRequest(
        {
            "Efield": electric_field,
            "cell": (np.diag([volume, 1.0, 1.0]), np.eye(3)),
            "result": (
                0.0,
                np.zeros(3),
                np.zeros((3, 3)),
                {
                    "dipole": dipole,
                    "BEC": np.zeros((1, 3, 3)),
                    "piezoelectric": np.zeros((3, 3, 3)),
                },
            ),
        }
    )

    dielectric.fixed_E(request)

    np.testing.assert_allclose(
        dielectric.get_electric_displacement(),
        [2.0 + 4.0 * np.pi * dipole[0] / volume, 0.0, 0.0],
    )


def _completed_request(extras):
    request = ForceRequest(
        {
            "id": 0,
            "status": "Done",
            "Efield": [0.0, 0.0, 0.1],
            "cell": (np.eye(3), np.eye(3)),
            "result": (0.0, np.zeros(3), np.zeros((3, 3)), extras),
        }
    )
    return request


def test_client_field_requires_an_applied_field_acknowledgement():
    dielectric = _client_dielectric()
    request = _completed_request({})
    dielectric.forcefield.requests.append(request)

    with pytest.raises(ValueError, match="applied_fields"):
        dielectric.post_process(request)


def test_client_field_warns_when_diagnostics_are_missing(capsys):
    dielectric = _client_dielectric()
    request = _completed_request({"applied_fields": {"electric_field": [0, 0, 0.1]}})
    dielectric.forcefield.requests.append(request)

    assert dielectric.post_process(request) is request
    assert "Please provide this value so that i-PI can verify" in (
        capsys.readouterr().out
    )


def test_client_field_rejects_missing_field_feedback():
    dielectric = _client_dielectric()
    request = _completed_request({"dipole": [0.0, 0.0, 0.0]})
    dielectric.forcefield.requests.append(request)

    with pytest.raises(ValueError, match="applied_fields"):
        dielectric.post_process(request)


def test_client_field_accepts_driver_acknowledgement():
    dielectric = _client_dielectric()
    request = _completed_request(
        {
            "dipole": [0.0, 0.0, 0.0],
            "applied_fields": {
                "electric_field": [0.0, 0.0, 0.1],
                "displacement_field": [0.0, 0.0, 0.1],
                "effective_electric_field": [0.0, 0.0, 0.1],
            },
        }
    )
    dielectric.forcefield.requests.append(request)

    assert dielectric.post_process(request) is request


def test_client_displacement_caches_field_quantity_for_properties():
    """Client-side fixed D exposes D - 4 pi mu / Omega as electric_field."""
    dielectric = _client_dielectric()
    volume = 2.0
    dipole = np.array([0.1, 0.0, 0.0])
    request = ForceRequest(
        {
            "id": 0,
            "status": "Done",
            "Dfield": [4.0 * np.pi * dipole[0] / volume + 2.0, 0.0, 0.0],
            "cell": (np.diag([volume, 1.0, 1.0]), np.eye(3)),
            "result": (
                0.0,
                np.zeros(3),
                np.zeros((3, 3)),
                {
                    "dipole": dipole,
                    "epsilon_infinity": np.eye(3),
                    "applied_fields": {
                        "electric_field": [2.0, 0.0, 0.0],
                        "displacement_field": [
                            4.0 * np.pi * dipole[0] / volume + 2.0,
                            0.0,
                            0.0,
                        ],
                        "effective_electric_field": [2.0, 0.0, 0.0],
                    },
                },
            ),
        }
    )
    dielectric.forcefield.requests.append(request)

    assert dielectric.post_process(request) is request
    np.testing.assert_array_equal(dielectric.get_electric_field(), [2.0, 0.0, 0.0])
