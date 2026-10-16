"""Normal-mode transformations preserve geometry and ring-polymer spring energy."""

import numpy as np
import pytest

from ipi.utils.nmtransform import (
    mk_nm_matrix,
    mk_o_nm_matrix,
    nm_eva,
    o_nm_eva,
    nm_fft,
    nm_trans,
    nm_rescale,
)


@pytest.mark.parametrize("nbeads", [1, 2, 3, 4, 7, 8])
@pytest.mark.parametrize("opened", [False, True])
def test_modes_diagonalize_the_bead_spring_energy(nbeads, opened):
    """The spring matrix is independently constructed from neighbouring beads."""
    transform = mk_o_nm_matrix(nbeads) if opened else mk_nm_matrix(nbeads)
    frequencies = o_nm_eva(nbeads) if opened else nm_eva(nbeads)
    spring = np.zeros((nbeads, nbeads))
    for bead in range(nbeads - 1 if opened else nbeads):
        difference = np.zeros(nbeads)
        difference[bead] += 1
        difference[(bead + 1) % nbeads] -= 1
        spring += np.outer(difference, difference)
    np.testing.assert_allclose(transform @ transform.T, np.eye(nbeads), atol=2e-14)
    np.testing.assert_allclose(
        transform @ spring @ transform.T, np.diag(frequencies**2), atol=2e-14
    )
    np.testing.assert_allclose(transform[0], np.full(nbeads, 1 / np.sqrt(nbeads)))


@pytest.mark.parametrize("nbeads", [1, 2, 3, 4, 7, 8])
@pytest.mark.parametrize("open_paths", [[], [0], [1]])
def test_fft_and_matrix_transforms_preserve_positions(nbeads, open_paths):
    """Odd/even bead counts and mixed open/closed paths use the same convention."""
    positions = np.random.default_rng(92).normal(size=(nbeads, 6))
    matrix = nm_trans(nbeads, open_paths=open_paths)
    fft = nm_fft(nbeads, 2, open_paths=open_paths)
    modes = matrix.b2nm(positions)
    np.testing.assert_allclose(fft.b2nm(positions), modes, atol=2e-14)
    np.testing.assert_allclose(fft.nm2b(modes), positions, atol=2e-14)
    np.testing.assert_allclose(matrix.nm2b(modes), positions, atol=2e-14)
    np.testing.assert_allclose(np.sum(modes**2), np.sum(positions**2), atol=2e-14)


@pytest.mark.parametrize("nsmall,nlarge", [(1, 4), (3, 8), (4, 7), (4, 4)])
@pytest.mark.parametrize("open_paths", [[], [0, 1]])
def test_mixed_path_expansion_preserves_each_atom(nsmall, nlarge, open_paths):
    """Expansion followed by contraction must recover both types of path."""
    positions = np.random.default_rng(17).normal(size=(nsmall, 6))
    rescale = nm_rescale(nsmall, nlarge, open_paths=open_paths)
    expanded = rescale.b1tob2(positions)
    np.testing.assert_allclose(
        expanded.mean(axis=0), positions.mean(axis=0), atol=2e-14
    )
    np.testing.assert_allclose(rescale.b2tob1(expanded), positions, atol=2e-14)
    assert not np.shares_memory(positions, expanded)


@pytest.mark.parametrize("open_paths", [[0], [1]])
@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="nm_rescale.b2tob1 applies the open-path map to closed-path atoms",
)
def test_expansion_preserves_mixed_open_and_closed_paths(open_paths):
    """An open path must not change the transformation of another atom's closed path."""
    positions = np.random.default_rng(17).normal(size=(3, 6))
    rescale = nm_rescale(3, 8, open_paths=open_paths)
    np.testing.assert_allclose(
        rescale.b2tob1(rescale.b1tob2(positions)), positions, atol=2e-14
    )
