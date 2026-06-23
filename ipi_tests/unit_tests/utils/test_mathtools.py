#!/usr/bin/env python3

import numpy as np
import numpy.testing as npt

from ipi.utils.mathtools import sinch


def test_sinch_accepts_small_complex_values():
    values = np.array([1.0 + 0.0j, 1.0e-8 + 1.0e-8j])

    expected = np.array(
        [
            np.sinh(values[0]) / values[0],
            1
            + (values[1] ** 2 / 6.0)
            * (1 + (values[1] ** 2 / 20.0) * (1 + (values[1] ** 2 / 42.0))),
        ]
    )

    npt.assert_allclose(sinch(values), expected)


def test_sinch_reconstruction_of_real_matrix_stays_real():
    v = np.array([[0.0, -1.0], [1.0, 0.0]])
    eigvals, eigvecs = np.linalg.eig(v)
    ieigvecs = np.linalg.inv(eigvecs)

    reconstructed = np.real_if_close(
        np.dot(eigvecs, np.dot(np.diag(sinch(eigvals)), ieigvecs))
    )

    npt.assert_allclose(reconstructed, np.sin(1.0) * np.eye(2))
    assert np.isrealobj(reconstructed)
