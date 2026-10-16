"""Numerical Hessians and rigid-motion removal have exact quadratic references."""

import numpy as np
import pytest
from ipi.utils.hesstools import clean_hessian, get_dynmat, get_hessian


@pytest.mark.parametrize("resume", [False, True])
@pytest.mark.parametrize("friction", [False, True])
def test_finite_difference_hessian_and_friction_derivative(
    tmp_path, monkeypatch, resume, friction
):
    monkeypatch.chdir(tmp_path)
    hessian = np.array([[2.0, 0.2, 0.1], [0.2, 3.0, -0.3], [0.1, -0.3, 4.0]])
    q = np.arange(6).reshape(2, 3) / 10

    class Quadratic:
        def __call__(self, x, new_disc=False):
            self.q = x.copy()
            self.eta = np.array([np.diag(1 + v**2) for v in x])
            return np.sum(x * (x @ hessian)) / 2, x @ hessian

    mapper = Quadratic()
    if resume:
        partial = np.zeros((3, 6))
        partial[0] = np.tile(hessian[0], 2)
        np.savetxt("hessian_0.tmp", partial)
        if friction:
            eta_partial = np.zeros((2, 3, 3, 3))
            eta_partial[:, 0, 0, 0] = 2 * q[:, 0]
            np.savetxt("hessianEta_0.tmp", eta_partial.ravel())
    result = get_hessian(mapper, q, 1, nbeads=2, friction=friction)
    if friction:
        result, eta_derivative = result
        expected = np.zeros((2, 3, 3, 3))
        for i in range(3):
            expected[:, i, i, i] = 2 * q[:, i]
        np.testing.assert_allclose(eta_derivative, expected, atol=1e-11)
    np.testing.assert_allclose(result, np.tile(hessian, (1, 2)), atol=1e-11)
    np.testing.assert_array_equal(mapper.q, q)
    assert not list(tmp_path.glob("*.tmp"))


@pytest.mark.parametrize("asr,removed", [("none", 0), ("crystal", 3), ("poly", 6)])
def test_rigid_motion_projection_preserves_remaining_unit_frequencies(asr, removed):
    q = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
    masses = np.array([1.0, 2.0, 3.0, 4.0])
    m3 = np.repeat(masses, 3)[None, :]
    hessian = np.diag(m3.ravel())
    eigenvalues, eigenvectors, inertia = clean_hessian(
        hessian, q, 4, 1, masses, m3, asr, mofi=True
    )
    assert eigenvalues.size == 12 - removed
    np.testing.assert_allclose(eigenvalues, 1, atol=1e-12)
    np.testing.assert_allclose(
        eigenvectors.T @ eigenvectors, np.eye(12 - removed), atol=1e-12
    )
    if removed:
        translations = np.kron(np.sqrt(masses)[:, None], np.eye(3))
        np.testing.assert_allclose(eigenvectors.T @ translations, 0, atol=1e-12)
    if asr == "poly":
        centered = q - np.average(q, axis=0, weights=masses)
        tensor = sum(
            m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
            for m, r in zip(masses, centered)
        )
        assert inertia == pytest.approx(np.linalg.det(tensor))
    np.testing.assert_allclose(get_dynmat(hessian, m3), np.eye(12))
