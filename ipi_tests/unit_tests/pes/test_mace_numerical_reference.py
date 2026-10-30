"""Real small MACE models: batching, force derivatives, and rotational covariance."""

import copy
import json

import numpy as np
import pytest

pytest.importorskip("mace")
torch = pytest.importorskip("torch")
from ase import Atoms
from ase.io import read, write
from e3nn import o3
from mace.modules import (
    MACE,
    RealAgnosticInteractionBlock,
    RealAgnosticResidualInteractionBlock,
)

from ipi.pes._mace import BatchedMACE, MACE_driver, run_cli
from ipi.utils.units import unit_to_internal


@pytest.fixture
def model():
    previous = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)
    try:
        with torch.random.fork_rng():
            torch.manual_seed(31)
            yield MACE(
                r_max=3.0,
                num_bessel=3,
                num_polynomial_cutoff=3,
                max_ell=1,
                interaction_cls=RealAgnosticResidualInteractionBlock,
                interaction_cls_first=RealAgnosticInteractionBlock,
                num_interactions=2,
                num_elements=1,
                hidden_irreps=o3.Irreps("4x0e + 4x1o"),
                MLP_irreps=o3.Irreps("4x0e"),
                atomic_energies=np.array([0.0]),
                avg_num_neighbors=2.0,
                atomic_numbers=[1],
                correlation=2,
                gate=torch.nn.functional.silu,
            )
    finally:
        torch.set_default_dtype(previous)


@pytest.mark.parametrize("batch_size", [1, 3])
def test_mace_batch_matches_single_and_energy_derivative(model, batch_size):
    calculator = BatchedMACE(
        models=[model],
        device="cpu",
        default_dtype="float64",
        instructions={"batch_size": batch_size, "ignore": ["contributions"]},
    )
    structures = [
        Atoms(
            "H" * n,
            positions=np.arange(n * 3).reshape(n, 3) * 0.25,
            cell=np.eye(3) * 8,
            pbc=True,
        )
        for n in [2, 3, 2]
    ]
    results = calculator.compute_batched(structures)
    for atoms, result in zip(structures, results):
        single = calculator.compute_batched([atoms])[0]
        for key in ["energy", "forces", "stress"]:
            np.testing.assert_allclose(result[key], single[key], atol=1e-10)
        np.testing.assert_allclose(result["forces"].sum(axis=0), 0, atol=1e-10)
    atoms = structures[0]
    delta = 1e-5
    plus, minus = atoms.copy(), atoms.copy()
    plus.positions[0, 0] += delta
    minus.positions[0, 0] -= delta
    perturbed = calculator.compute_batched([plus, minus])
    assert results[0]["forces"][0, 0] == pytest.approx(
        -(perturbed[0]["energy"] - perturbed[1]["energy"]) / (2 * delta), abs=1e-7
    )
    # A proper rotation changes force components but not energy.
    rotation = np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, 1]])
    rotated = atoms.copy()
    rotated.positions[:] = atoms.positions @ rotation.T
    rotated.cell[:] = atoms.cell.array @ rotation.T
    transformed = calculator.compute_batched([rotated])[0]
    assert transformed["energy"] == pytest.approx(results[0]["energy"], abs=1e-10)
    np.testing.assert_allclose(
        transformed["forces"], results[0]["forces"] @ rotation.T, atol=1e-10
    )


def test_identical_mace_committee_matches_single_model(model):
    atoms = Atoms(
        "HH", positions=[[0, 0, 0], [0.8, 0.2, 0]], cell=np.eye(3) * 8, pbc=True
    )
    single = BatchedMACE(
        models=[copy.deepcopy(model)],
        device="cpu",
        default_dtype="float64",
        instructions={"ignore": ["contributions"]},
    ).compute_batched([atoms])[0]
    committee = BatchedMACE(
        models=[model, copy.deepcopy(model)],
        device="cpu",
        default_dtype="float64",
        instructions={"ignore": ["contributions"]},
    ).compute_batched([atoms])[0]
    for key in ["energy", "forces", "stress"]:
        np.testing.assert_allclose(committee[key], single[key], atol=1e-12)


def test_mace_cli_exports_model_results_with_requested_prefix(
    model, tmp_path, monkeypatch
):

    checkpoint = tmp_path / "model.pt"
    torch.save(model, checkpoint)
    atoms = Atoms(
        "HH", positions=[[0.1, 0, 0], [0.8, 0.2, 0]], cell=np.eye(3) * 8, pbc=True
    )
    source = tmp_path / "input.xyz"
    output = tmp_path / "output.xyz"
    write(source, [atoms, atoms.copy()])
    settings = tmp_path / "settings.json"
    settings.write_text(
        json.dumps(
            {"default_dtype": "float64", "instructions": {"ignore": ["contributions"]}}
        )
    )
    shapes = tmp_path / "shapes.json"
    shapes.write_text(json.dumps({"forces": ["natoms", 3]}))
    monkeypatch.setattr(
        "sys.argv",
        [
            "mace",
            "-m",
            str(checkpoint),
            "-i",
            str(source),
            "-o",
            str(output),
            "-mk",
            str(settings),
            "--ase-like-properties",
            str(shapes),
            "-p",
            "test_",
        ],
    )
    run_cli()
    frames = read(output, index=":")
    assert len(frames) == 2
    np.testing.assert_array_equal(frames[0].positions, atoms.positions)
    np.testing.assert_allclose(
        frames[0].arrays["test_forces"], frames[1].arrays["test_forces"]
    )
    driver = MACE_driver(str(source), str(checkpoint), mace_kwargs=str(settings))
    cell = unit_to_internal("length", "angstrom", atoms.cell.array.T)
    q = unit_to_internal("length", "angstrom", atoms.positions)
    for count in [1, 3, 2]:
        results = driver.compute(
            [cell.copy() for _ in range(count)], [q.copy() for _ in range(count)]
        )
        if count == 1:
            results = [results]
        assert len(results) == count
        for energy, force, virial, _ in results:
            assert energy == pytest.approx(
                unit_to_internal(
                    "energy", "electronvolt", frames[0].info["test_energy"]
                ),
                abs=1e-10,
            )
            np.testing.assert_allclose(
                force,
                unit_to_internal("force", "ev/ang", frames[0].arrays["test_forces"]),
                atol=1e-9,
            )
            np.testing.assert_allclose(
                virial,
                unit_to_internal(
                    "energy",
                    "electronvolt",
                    -frames[0].info["test_stress"] * atoms.get_volume(),
                ),
                atol=1e-9,
            )
    with pytest.raises(ValueError, match="volume"):
        driver.template2atoms([np.eye(3) * 1e5], [atoms.positions])


@pytest.mark.parametrize(
    "shapes,error",
    [
        (None, TypeError),
        ({1: []}, TypeError),
        ({"x": 3}, TypeError),
        ({"x": [-1]}, ValueError),
        ({"x": [3, "natoms"]}, ValueError),
    ],
)
def test_mace_rejects_invalid_output_shapes(shapes, error):
    with pytest.raises(error):
        BatchedMACE._normalize_ase_like_properties(shapes)
