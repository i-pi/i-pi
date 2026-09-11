"""Tests for containers used to organize potential-energy results."""

import numpy as np
import pytest

from ipi.pes.tools import ModelResults, StructureResults


def test_structure_results_stores_scalar_and_atomic_properties():
    """Checks scalar storage and expansion of per-atom result shapes."""

    results = StructureResults(natoms=2, shapes={"energy": (), "forces": ("natoms", 3)})
    forces = np.arange(6, dtype=float)

    results.store("energy", 1.5)
    results.store("forces", forces)

    assert results.shapes == {"energy": (), "forces": (2, 3)}
    assert results["energy"] == 1.5
    assert np.array_equal(results["forces"], forces.reshape(2, 3))


def test_structure_results_rejects_unknown_or_invalid_properties():
    """Checks that invalid result names and shapes are rejected."""

    results = StructureResults(natoms=2, shapes={"forces": ("natoms", 3)})

    with pytest.raises(KeyError, match="Unknown property"):
        results.store("energy", 1.5)
    with pytest.raises(ValueError, match="wrong shape"):
        results.store("forces", np.zeros((3, 3)))


def test_structure_results_dictionary_array_copy_is_optional():
    """Checks the documented copied and shared array views."""

    results = StructureResults(natoms=2, shapes={"forces": ("natoms", 3)})
    results.store("forces", np.arange(6, dtype=float))

    copied = results.as_dict()
    copied["forces"][0, 0] = -1.0
    assert results["forces"][0, 0] == 0.0

    shared = results.as_dict(copy_arrays=False)
    shared["forces"][0, 0] = -2.0
    assert results["forces"][0, 0] == -2.0


def test_model_results_splits_atomic_results_between_structures():
    """Checks batched outputs for structures with different atom counts."""

    results = ModelResults({"energy": (), "forces": ("natoms", 3)})
    energies = np.array([1.0, 2.0])
    forces = np.arange(9, dtype=float).reshape(3, 3)

    results.store([1, 2], {"energy": energies, "forces": forces})

    assert len(results) == 2
    assert results.natoms() == [1, 2]
    assert results[0]["energy"] == energies[0]
    assert results[1]["energy"] == energies[1]
    assert np.array_equal(results[0]["forces"], forces[:1])
    assert np.array_equal(results[1]["forces"], forces[1:])


def test_model_results_mean_averages_compatible_models():
    """Checks averaging of scalar and per-atom model predictions."""

    shapes = {"energy": (), "forces": ("natoms", 3)}
    first = ModelResults(shapes)
    second = ModelResults(shapes)
    first.store(
        [1, 2],
        {"energy": [1.0, 3.0], "forces": np.arange(9, dtype=float).reshape(3, 3)},
    )
    second.store(
        [1, 2],
        {
            "energy": [3.0, 5.0],
            "forces": np.arange(9, 18, dtype=float).reshape(3, 3),
        },
    )

    averaged = ModelResults.mean([first, second])

    assert averaged.natoms() == [1, 2]
    assert averaged[0]["energy"] == 2.0
    assert averaged[1]["energy"] == 4.0
    assert np.array_equal(
        averaged[0]["forces"], (first[0]["forces"] + second[0]["forces"]) / 2.0
    )
    assert np.array_equal(
        averaged[1]["forces"], (first[1]["forces"] + second[1]["forces"]) / 2.0
    )


def test_model_results_mean_rejects_empty_or_incompatible_models():
    """Checks that a model average requires matching structures."""

    with pytest.raises(ValueError, match="empty"):
        ModelResults.mean([])

    shapes = {"energy": ()}
    first = ModelResults(shapes)
    incompatible = ModelResults(shapes)
    first.store([1], {"energy": [1.0]})
    incompatible.store([2], {"energy": [2.0]})

    with pytest.raises(ValueError, match="same natoms and shapes"):
        ModelResults.mean([first, incompatible])
