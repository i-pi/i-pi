"""Tests for FFDielectric field application checks."""

import numpy as np
import pytest

from ipi.engine.forcefields import FFDielectric, ForceField, ForceRequest


def _client_dielectric():
    """Create the minimum FFDielectric object needed for result validation."""
    wrapped = ForceField()
    return FFDielectric(
        name="dielectric",
        where="client",
        dipole={"family": "dipole", "units": "atomic_unit", "key": "dipole"},
        bec={"family": "born_charge", "units": "atomic_unit", "key": "BEC"},
        piezo={"family": "electric-polarization", "units": "atomic_unit", "key": "piezoelectric"},
        electric_fields=[],
        electric_displacements=[],
        forcefield=wrapped,
    )


def _completed_request(extras):
    request = ForceRequest(
        {
            "id": 0,
            "status": "Done",
            "Efield": [0.0, 0.0, 0.1],
            "result": (0.0, np.zeros(3), np.zeros((3, 3)), extras),
        }
    )
    return request


def test_client_field_requires_driver_acknowledgement():
    dielectric = _client_dielectric()
    request = _completed_request({})
    dielectric.forcefield.requests.append(request)

    with pytest.raises(ValueError, match="applied_fields"):
        dielectric.post_process(request)


def test_client_field_accepts_driver_acknowledgement():
    dielectric = _client_dielectric()
    request = _completed_request({"applied_fields": ["electric_field"]})
    dielectric.forcefield.requests.append(request)

    assert dielectric.post_process(request) is request
