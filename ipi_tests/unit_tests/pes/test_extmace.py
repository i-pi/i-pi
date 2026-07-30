"""Tests for the FFDielectric additions to the MACE driver."""

import json

import pytest

pytest.importorskip("mace")

from ipi.pes._mace import MACE_driver
from ipi.pes.extmace import Extended_MACE_driver


def test_extmace_acknowledges_applied_electric_field(monkeypatch):
    """An Efield request is acknowledged in the returned extras JSON."""
    result = (1.0, "forces", "virial", '{"dipole": [0.0, 0.0, 0.0]}')
    monkeypatch.setattr(MACE_driver, "post_process", lambda *args: result)

    driver = object.__new__(Extended_MACE_driver)
    driver.extra = {"Efield": [0.0, 0.0, 0.1]}

    _, _, _, extras = driver.post_process({}, None)

    assert json.loads(extras)["applied_fields"] == ["electric_field"]
