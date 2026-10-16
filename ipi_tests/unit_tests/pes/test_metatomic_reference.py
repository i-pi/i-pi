"""Exported quadratic pair model checks metatomic forces, virials, and batching."""

from typing import Dict, List, Optional

import numpy as np
import pytest

torch = pytest.importorskip("torch")
mta = pytest.importorskip("metatomic.torch")
pytest.importorskip("vesin.metatomic")
from metatensor.torch import Labels, TensorBlock, TensorMap
from metatomic.torch import System, ModelOutput
from ase import Atoms
from ase.io import write
from ipi.pes.metatomic import MetatomicDriver
from ipi.utils.units import unit_to_internal, unit_to_user


class PairQuadratic(torch.nn.Module):
    """Energy in eV is half the squared pair separation in angstroms."""

    def forward(
        self,
        systems: List[System],
        outputs: Dict[str, ModelOutput],
        selected_atoms: Optional[Labels] = None,
    ) -> Dict[str, TensorMap]:
        values = torch.stack(
            [
                0.5 * torch.sum((system.positions[1] - system.positions[0]) ** 2)
                for system in systems
            ]
        ).reshape(-1, 1)
        block = TensorBlock(
            values=values,
            samples=Labels(
                ["system"], torch.arange(len(systems), dtype=torch.int32).reshape(-1, 1)
            ),
            components=torch.jit.annotate(List[Labels], []),
            properties=Labels(["energy"], torch.zeros((1, 1), dtype=torch.int32)),
        )
        result = {
            "energy": TensorMap(
                Labels(["_"], torch.zeros((1, 1), dtype=torch.int32)), [block]
            )
        }
        if "energy_ensemble" in outputs:
            ensemble_block = TensorBlock(
                values=values * torch.tensor([[0.8, 1.2]], dtype=values.dtype),
                samples=Labels(
                    ["system"],
                    torch.arange(len(systems), dtype=torch.int32).reshape(-1, 1),
                ),
                components=torch.jit.annotate(List[Labels], []),
                properties=Labels(
                    ["energy"], torch.tensor([[0], [1]], dtype=torch.int32)
                ),
            )
            result["energy_ensemble"] = TensorMap(
                Labels(["_"], torch.zeros((1, 1), dtype=torch.int32)), [ensemble_block]
            )
        return result


@pytest.fixture
def model_files(tmp_path):
    model = mta.AtomisticModel(
        PairQuadratic().eval(),
        mta.ModelMetadata(name="quadratic test"),
        mta.ModelCapabilities(
            outputs={
                "energy": mta.ModelOutput(quantity="energy", unit="eV", per_atom=False),
                "energy_ensemble": mta.ModelOutput(
                    quantity="energy", unit="eV", per_atom=False
                ),
            },
            atomic_types=[1],
            interaction_range=0.0,
            length_unit="angstrom",
            supported_devices=["cpu"],
            dtype="float64",
        ),
    )
    path = tmp_path / "pair.pt"
    model.save(str(path))
    template = tmp_path / "template.xyz"
    write(template, Atoms("HH", positions=[[0, 0, 0], [1, 0, 0]]))
    return str(template), str(path)


@pytest.mark.parametrize("devices", ["cpu", "cpu;cpu"])
def test_metatomic_pair_force_virial_and_batch_order(model_files, devices):
    driver = MetatomicDriver(*model_files, device=devices, check_consistency=True)
    positions = [
        np.array([[0.1, 0.2, 0.3], [0.4, -0.1, 0.2]]) * scale
        for scale in [1.0, 2.0, 0.5]
    ]
    cells = [np.eye(3) * 10 for _ in positions]
    try:
        batch = driver.compute(cells, positions)
        for q, cell, result in zip(positions, cells, batch):
            energy, force, virial, extra = result
            dr = (q[1] - q[0]) * unit_to_user("length", "angstrom", 1)
            assert energy == pytest.approx(
                unit_to_internal("energy", "electronvolt", 0.5 * np.dot(dr, dr))
            )
            np.testing.assert_allclose(
                force,
                unit_to_internal("force", "ev/ang", np.array([dr, -dr])),
                atol=1e-12,
            )
            np.testing.assert_allclose(
                virial,
                unit_to_internal("energy", "electronvolt", -np.outer(dr, dr)),
                atol=1e-12,
            )
            assert extra == "{}"
            single = driver.compute(cell, q)
            for k in range(3):
                np.testing.assert_allclose(single[k], result[k], atol=1e-12)
        with pytest.raises(ValueError, match="Both position and cell"):
            driver.compute(cells, positions[:1])
    finally:
        if driver.executor is not None:
            driver.executor.shutdown()


@pytest.mark.parametrize(
    "options",
    [
        {"force_virial_ensemble": True},
        {
            "force_virial_ensemble": True,
            "energy_ensemble": True,
            "non_conservative": True,
        },
    ],
)
def test_metatomic_rejects_conflicting_ensemble_options(options):
    with pytest.raises(ValueError):
        MetatomicDriver("unused.xyz", "unused.pt", **options)


def test_metatomic_committee_derivatives_match_each_quadratic(model_files):
    import json

    driver = MetatomicDriver(
        *model_files,
        energy_ensemble=True,
        force_virial_ensemble=True,
        check_consistency=True
    )
    cell = np.eye(3) * 10
    q = np.array([[0.1, 0.2, 0.3], [0.5, -0.2, 0.1]])
    energy, force, virial, raw = driver.compute(cell, q)
    extra = json.loads(raw)
    weights = np.array([0.8, 1.2])
    np.testing.assert_allclose(extra["committee_pot"], weights * energy, atol=1e-12)
    np.testing.assert_allclose(
        np.array(extra["committee_force"]).reshape(2, 2, 3),
        weights[:, None, None] * force,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        np.array(extra["committee_virial"]).reshape(2, 3, 3),
        weights[:, None, None] * virial,
        atol=1e-12,
    )
