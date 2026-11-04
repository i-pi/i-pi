"""Tests for server-side fixed-electric-field FFDielectric support."""

import numpy as np
import pytest

from ipi.engine.forcefields import FFDielectric, ForceField, ForceRequest


def _server_dielectric(**kwargs):
    """Create the smallest server-side dielectric wrapper for unit tests."""
    return FFDielectric(
        name="dielectric",
        dipole={"family": "electric-dipole", "units": "atomic_unit", "key": "dipole"},
        bec={"family": "charge", "units": "e", "key": "BEC"},
        piezo={
            "family": "electric-polarization",
            "units": "atomic_unit",
            "key": "piezoelectric",
        },
        electric_fields=[],
        forcefield=ForceField(),
        **kwargs,
    )


def test_dielectric_is_server_only_for_now():
    """Client-side coupling belongs to the follow-up PR."""
    with pytest.raises(ValueError, match="where='server'"):
        _server_dielectric(where="client")


def test_dielectric_rejects_wrapped_forcefield():
    """The wrapper must receive unwrapped coordinates."""
    with pytest.raises(ValueError, match="does not support pbc='True'"):
        FFDielectric(
            name="dielectric",
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
            forcefield=ForceField(dopbc=True),
        )


def test_dielectric_applies_fixed_e_energy_force_and_virial():
    """Fixed E uses the driver dipole, BEC, and proper piezoelectric tensor."""
    dielectric = _server_dielectric()
    field = np.array([2.0, 0.0, 0.0])
    piezoelectric = np.zeros((3, 3, 3))
    piezoelectric[0, 0, 0] = 4.0
    request = ForceRequest(
        {
            "electric_field": field,
            "cell": (np.diag([2.0, 1.0, 1.0]), np.eye(3)),
            "result": (
                10.0,
                np.zeros(3),
                np.zeros((3, 3)),
                {
                    "dipole": np.array([0.5, 0.0, 0.0]),
                    "BEC": np.array(
                        [[[3.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]
                    ),
                    "piezoelectric": piezoelectric,
                },
            ),
        }
    )

    energy, forces, virial, _ = dielectric.fixed_E(request)

    assert energy == pytest.approx(9.0)
    np.testing.assert_allclose(forces, [6.0, 0.0, 0.0])
    np.testing.assert_allclose(
        virial, [[16.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    )


def test_dielectric_expands_voigt_piezoelectric_tensor():
    """Voigt piezoelectric data expands to symmetric Cartesian strain axes."""
    voigt = np.arange(18).reshape(3, 6)
    tensor = FFDielectric._coerce_piezoelectric_tensor(voigt)
    np.testing.assert_array_equal(tensor[:, 1, 2], voigt[:, 3])
    np.testing.assert_array_equal(tensor[:, 0, 2], voigt[:, 4])
    np.testing.assert_array_equal(tensor[:, 0, 1], voigt[:, 5])
    np.testing.assert_array_equal(tensor, tensor.swapaxes(1, 2))


def test_dielectric_rejects_nonsymmetric_piezoelectric_tensor():
    """The proper piezoelectric tensor must be symmetric in strain axes."""
    tensor = np.zeros((3, 3, 3))
    tensor[0, 1, 2] = 1.0
    with pytest.raises(ValueError, match="symmetric in its two strain indices"):
        FFDielectric._coerce_piezoelectric_tensor(tensor)
