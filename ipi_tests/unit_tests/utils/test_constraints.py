"""Constraint gradients and projections checked against molecular geometry."""

import numpy as np
import pytest

from ipi.engine.beads import Beads
from ipi.engine.motion.constrained_dynamics import ConstraintSolver
from ipi.utils.constrtools import (
    AngleConstraint,
    ConstraintList,
    EckartConstraint,
    RigidBondConstraint,
)


@pytest.fixture
def molecule():
    """A noncollinear molecule with unequal masses and nonzero momenta."""
    beads = Beads(3, 1)
    beads.m[:] = [2.0, 3.0, 5.0]
    beads.q[:] = [[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]]
    beads.p[:] = [[0.1, 0.2, 0.3, -0.3, 0.4, 0.2, 0.5, -0.2, 0.1]]
    return beads


def make_constraint(kind, molecule, infer=False):
    """Construct the geometry constraints without evaluating their implementation."""
    if kind == "bonds":
        constraint = RigidBondConstraint(
            np.array([[0, 1], [0, 2]]), np.array([] if infer else [1.0, 1.0])
        )
    elif kind == "angle":
        constraint = AngleConstraint(
            np.array([[0, 1, 2]]), np.array([] if infer else [np.pi / 2])
        )
    elif kind == "eckart":
        constraint = EckartConstraint(
            np.array([0, 1, 2]), np.array([]) if infer else molecule.q[0].copy()
        )
    else:
        constraint = ConstraintList(
            [make_constraint("bonds", molecule), make_constraint("angle", molecule)]
        )
    return constraint


@pytest.mark.parametrize("kind", ["bonds", "angle", "eckart", "combined"])
@pytest.mark.parametrize("infer", [False, True])
def test_constraint_jacobian_matches_finite_displacements(molecule, kind, infer):
    """The Jacobian must predict changes in every constraint, including combined lists."""
    constraint = make_constraint(kind, molecule, infer)
    constraint.bind(molecule)
    reference = molecule.q[0, constraint.i3_unique].copy()
    constraint.q = reference
    constraint.qprev = reference
    np.testing.assert_allclose(constraint.g, 0, atol=1e-14)
    analytic = constraint.Dg.copy()
    numerical = np.zeros_like(analytic)
    delta = 1e-6
    for component in range(reference.size):
        constraint.q = reference.copy()
        constraint.q[component] += delta
        forward = constraint.g.copy()
        constraint.q[component] -= 2 * delta
        backward = constraint.g.copy()
        numerical[:, component] = (forward - backward) / (2 * delta)
    np.testing.assert_allclose(analytic, numerical, rtol=1e-8, atol=1e-9)
    gram = constraint.Gram
    np.testing.assert_allclose(gram, gram.T, atol=1e-14)
    assert np.linalg.eigvalsh(gram).min() > 0


@pytest.mark.parametrize("kind", ["bonds", "angle", "eckart", "combined"])
def test_constraint_solver_projects_positions_and_velocities(molecule, kind):
    """After projection, positions satisfy the constraints and velocities preserve them."""
    constraint = make_constraint(kind, molecule)
    constraint.bind(molecule)
    solver = ConstraintSolver([constraint], dt=0.1, tolerance=1e-12)
    solver.bind(molecule)
    molecule.q[:] += np.random.default_rng(27).normal(scale=0.01, size=(1, 9))
    solver.proj_manifold()
    constraint.q = molecule.q[0, constraint.i3_unique].copy()
    np.testing.assert_allclose(constraint.g, 0, atol=1e-12)
    solver.proj_cotangent()
    velocity = molecule.p[0, constraint.i3_unique] / constraint.m3
    np.testing.assert_allclose(constraint.Dg @ velocity, 0, atol=1e-12)
    previous = molecule.p.copy()
    solver.proj_cotangent()
    np.testing.assert_allclose(molecule.p, previous, atol=1e-12)


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="Constraint Gram matrix is not invalidated by changed masses",
)
def test_constraint_gram_updates_when_masses_change(molecule):
    """Doubling masses halves the inverse-mass-weighted constraint matrix."""
    constraint = make_constraint("bonds", molecule)
    constraint.bind(molecule)
    constraint.qprev = molecule.q[0].copy()
    original = constraint.Gram.copy()
    molecule.m[:] *= 2
    np.testing.assert_allclose(constraint.Gram, original / 2, atol=1e-14)


def test_constraint_solver_rejects_path_integrals():
    """The solver explicitly supports classical dynamics only."""
    with pytest.raises(ValueError, match="classical MD"):
        ConstraintSolver([]).bind(Beads(1, 2))
