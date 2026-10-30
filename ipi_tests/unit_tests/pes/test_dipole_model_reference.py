"""A charge-linear Torch model checks dipole/BEC loading, units, and serialization."""

import json

import numpy as np
import pytest

torch = pytest.importorskip("torch")
from ipi.pes.driverdipole import driverdipole_driver, get_class


class ChargeDipole(torch.nn.Module):
    """Minimal physical model satisfying the documented dipole-model interface."""

    def __init__(self):
        super().__init__()
        self.register_buffer("charges", torch.tensor([1.0, -1.0], dtype=torch.float64))

    def n_parameters(self):
        return 0

    def get(self, cell, pos):
        q = torch.as_tensor(pos, dtype=torch.float64)
        return (self.charges[:, None] * q).sum(dim=0), None

    def get_value_and_jac(self, cell, pos):
        q = torch.as_tensor(pos, dtype=torch.float64)
        jac = torch.autograd.functional.jacobian(
            lambda r: (self.charges[:, None] * r).sum(dim=0), q
        )
        return self.get(cell, pos)[0], jac.reshape(3, -1).T, None


@pytest.mark.parametrize("bec", [False, True])
@pytest.mark.parametrize("send", [False, True])
def test_dipole_model_checkpoint_and_born_charges(tmp_path, bec, send):
    instructions = tmp_path / "model.json"
    instructions.write_text(
        json.dumps(
            {
                "module": __name__,
                "class": "ChargeDipole",
                "kwargs": {},
                "chemical-symbols": ["H", "H"],
            }
        )
    )
    weights = tmp_path / "model.pth"
    torch.save(ChargeDipole().state_dict(), weights)
    options = tmp_path / "options.json"
    dipole_file = tmp_path / "dipole.dat"
    bec_file = tmp_path / "bec.dat"
    options.write_text(
        json.dumps(
            {
                "dipole": {"send": send, "file": str(dipole_file)},
                "BEC": {"compute": bec, "send": send, "file": str(bec_file)},
            }
        )
    )
    driver = driverdipole_driver(False, str(instructions), str(weights), str(options))
    q = np.array([[0.2, -0.1, 0.4], [-0.3, 0.5, 0.1]])
    energy, force, virial, extra = driver.compute_structure(np.eye(3) * 10, q)
    assert energy == 0
    np.testing.assert_array_equal(force, np.zeros_like(q))
    np.testing.assert_array_equal(virial, np.zeros((3, 3)))
    output = json.loads(extra)
    if send:
        np.testing.assert_allclose(output["dipole"], q[0] - q[1])
        if bec:
            np.testing.assert_allclose(
                output["BEC"], np.vstack([np.eye(3), -np.eye(3)])
            )
    else:
        assert output == {}
    np.testing.assert_allclose(np.loadtxt(dipole_file), q[0] - q[1])
    if bec:
        np.testing.assert_allclose(
            np.loadtxt(bec_file), np.vstack([np.eye(3), -np.eye(3)])
        )
    assert driver.model._symbols == ["H", "H"]
    assert driver.count == 1


@pytest.mark.parametrize(
    "module,cls", [("no_such_ipi_model", "Model"), ("numpy", "NoSuchModel")]
)
def test_dipole_model_reports_missing_module_or_class(module, cls):
    with pytest.raises(ValueError, match="not found"):
        get_class(module, cls)
