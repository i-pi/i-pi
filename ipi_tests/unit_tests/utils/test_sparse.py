"""Tests for i-PI's sparse matrix storage helpers."""

import numpy as np
import pytest

from ipi.utils.sparse import csc_matrix, csr_matrix, load_npz, save_npz


@pytest.mark.parametrize("matrix_type", [csr_matrix, csc_matrix])
def test_sparse_matrix_roundtrip_to_dense(matrix_type):
    """Checks that both sparse formats preserve values and density."""

    dense = np.array([[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]])
    sparse = matrix_type(dense)

    assert np.array_equal(sparse.toarray(), dense)
    assert sparse.density() == pytest.approx(np.count_nonzero(dense) / dense.size)


def test_sparse_matrix_arithmetic_matches_numpy():
    """Checks sparse arithmetic against the equivalent dense operations."""

    first_dense = np.array([[1.0, 0.0], [2.0, 3.0]])
    second_dense = np.array([[4.0, 5.0], [0.0, 6.0]])
    first = csr_matrix(first_dense)
    second = csc_matrix(second_dense)

    assert np.array_equal((first + second).toarray(), first_dense + second_dense)
    assert np.array_equal((first - second).toarray(), first_dense - second_dense)
    assert np.array_equal((first * second).toarray(), first_dense * second_dense)
    assert np.array_equal(first.dot(second).toarray(), first_dense.dot(second_dense))


@pytest.mark.parametrize("matrix_type", [csr_matrix, csc_matrix])
@pytest.mark.parametrize("compressed", [True, False])
def test_sparse_npz_roundtrip_preserves_format_and_values(
    tmp_path, matrix_type, compressed
):
    """Checks saving and loading both sparse formats."""

    dense = np.array([[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]])
    original = matrix_type(dense)
    path = tmp_path / "matrix.npz"

    save_npz(path, original, compressed=compressed)
    restored = load_npz(path)

    assert restored.kind == original.kind
    assert np.array_equal(restored.toarray(), dense)
