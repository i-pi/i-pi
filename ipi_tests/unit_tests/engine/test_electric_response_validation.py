"""Validate force-client dipole and effective-charge data at the input boundary."""

from types import SimpleNamespace
import numpy as np
import pytest
from ipi.engine.beads import Beads
from ipi.engine.motion.driven_dynamics import BEC, ElectricDipole
from ipi.utils.depend import depend_value
from ipi.utils.units import unit_to_internal


def response_objects(extras):
    # Only the external client's response is supplied; parsing uses real objects.
    ensemble = SimpleNamespace(beads=Beads(2, 2), forces=SimpleNamespace(extras=extras))
    bec = BEC(cbec=True)
    bec.bind(ensemble, "eda-nve", depend_value(name="threshold", value=1e-8))
    dipole = ElectricDipole()
    dipole.bind(ensemble)
    return bec, dipole


def test_client_born_charges_obey_neutral_translation_sum_rule():
    charges = np.vstack([np.eye(3), -np.eye(3)])
    bec, _ = response_objects({"BEC": [charges.tolist(), (2 * charges).tolist()]})
    np.testing.assert_allclose(bec.bec, [charges, 2 * charges])
    translation = np.tile([0.3, -0.2, 0.5], 2)
    np.testing.assert_allclose(translation @ bec.bec, 0, atol=1e-15)


@pytest.mark.parametrize(
    "data,match",
    [
        ({}, "not returned"),
        ({"BEC": [np.zeros((3, 3))] * 2}, "number of BEC"),
        ({"BEC": [np.zeros((6, 2))] * 2}, "wrong shape"),
        ({"BEC": [np.ones((6, 3))] * 2}, "charge conservation"),
    ],
)
def test_invalid_client_charges_are_rejected(data, match):
    bec, _ = response_objects(data)
    with pytest.raises(ValueError, match=match):
        _ = bec.bec


@pytest.mark.parametrize("index", [-1, 2])
def test_response_bead_selection_rejects_invalid_indices(index):
    bec, dipole = response_objects({})
    with pytest.raises(ValueError, match="bead"):
        bec._get_driver_BEC(index)
    with pytest.raises(ValueError, match="bead"):
        dipole._get_dipole(index)


def test_structured_dipoles_preserve_bead_order():
    values = [[0.1, -0.2, 0.3], [0.4, 0.5, -0.6]]
    _, dipole = response_objects({"dipole": values})
    np.testing.assert_allclose(dipole.dipole, values)


@pytest.mark.parametrize(
    "extras", [{}, {"dipole": [[1, 2], [3, 4]]}, {"raw": ["invalid", "invalid"]}]
)
def test_missing_or_malformed_dipoles_cannot_look_like_zero(extras):
    _, dipole = response_objects(extras)
    assert np.isnan(dipole.dipole).all()


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="Legacy dipole parser multiplies a string array by a float and returns NaN",
)
def test_legacy_dipole_text_converts_angstrom_to_atomic_units():
    _, dipole = response_objects(
        {"raw": ["Total dipole moment [eAng]: 1.0 -2.0 3.0"] * 2}
    )
    expected = np.tile([1, -2, 3], (2, 1)) * unit_to_internal("length", "angstrom", 1)
    np.testing.assert_allclose(dipole.dipole, expected)
