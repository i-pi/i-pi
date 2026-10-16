"""Tests ring polymer contraction."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import pytest
import numpy as np
from numpy.testing import assert_almost_equal

from ipi.utils import nmtransform

numbers_to_check = list(range(10, 56, 9))


qs_to_check = [
    np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]]),
    np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], [0.0, 0.1, 0.0, 1.0, 0.1, 0.0]]),
    np.array(
        [
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.1, 0.0, 1.0, 0.1, 0.0],
            [0.0, -0.1, 0.0, 1.0, -0.1, 0.0],
        ]
    ),
    np.array(
        [
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.1, 0.0, 1.0, 0.1, 0.0],
            [0.0, 0.2, 0.0, 1.0, 0.2, 0.0],
            [0.0, -0.1, 0.0, 1.0, -0.1, 0.0],
        ]
    ),
]


# Preserve the original geometries and bead counts, then add open-path cases.
CASES = (
    [
        pytest.param(q, n, None, id=f"original-p{len(q)}-to-{n}")
        for q in qs_to_check
        for n in numbers_to_check
    ]
    + [
        pytest.param(
            np.random.default_rng(17).normal(size=(small, 6)),
            large,
            opened,
            id=f"{label}-p{small}-to-{large}",
        )
        for small, large in [(1, 4), (3, 8), (4, 7), (4, 4)]
        for label, opened in [("closed", []), ("open-both", [0, 1])]
    ]
    + [
        pytest.param(
            np.random.default_rng(17).normal(size=(3, 6)),
            8,
            [atom],
            id=f"mixed-open{atom}-p3-to-8",
            marks=pytest.mark.xfail(
                strict=True,
                raises=AssertionError,
                reason="nm_rescale.b2tob1 applies the open-path map to closed-path atoms",
            ),
        )
        for atom in [0, 1]
    ]
)


@pytest.mark.parametrize("q,n,open_paths", CASES)
def test_contraction(q, n, open_paths):
    """Expansion preserves each path, its centroid, and independent storage."""
    rescale = nmtransform.nm_rescale(len(q), n, open_paths=open_paths)
    expanded = rescale.b1tob2(q)
    restored = rescale.b2tob1(expanded)
    reverse = nmtransform.nm_rescale(n, len(q), open_paths=open_paths)

    assert_almost_equal(q, restored)
    assert_almost_equal(restored, reverse.b1tob2(expanded))
    assert_almost_equal(
        nmtransform.mk_rs_matrix(n, 1) @ expanded,
        nmtransform.mk_rs_matrix(len(q), 1) @ q,
    )
    np.testing.assert_allclose(expanded.mean(axis=0), q.mean(axis=0), atol=2e-14)
    np.testing.assert_allclose(restored, q, atol=2e-14)
    assert not np.shares_memory(q, expanded)
