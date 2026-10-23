"""Reference formulas for path-integral tests, written with NumPy only.

SC weights: Kapil, Behler & Ceriotti, arXiv:1606.00920, Eqs. 7-8,
Appendix C. Scaled-coordinate derivatives: Yamamoto, physics/0505109,
Eqs. 18,22,33,34. Resizing: Fourier projection, Markland & Manolopoulos,
JCP 129,024105 (2008).
"""

import itertools

import numpy as np
from numpy.polynomial import Polynomial


def fourier_matrix(count):
    """Real orthonormal Fourier basis in i-PI's documented mode order."""
    beads = np.arange(count)
    rows = []
    for mode in range(count):
        if mode == 0:
            row = np.ones(count) / np.sqrt(count)
        elif count % 2 == 0 and mode == count // 2:
            row = (-1.0) ** beads / np.sqrt(count)
        elif mode <= count // 2:
            row = np.sqrt(2 / count) * np.cos(2 * np.pi * mode * beads / count)
        else:
            row = np.sqrt(2 / count) * np.sin(2 * np.pi * mode * beads / count)
        rows.append(row)
    return np.asarray(rows)


def projection(source, target):
    """Coordinate projection for odd smaller paths (no Nyquist ambiguity)."""
    small = min(source, target)
    if source == target:
        return np.eye(source)
    if small % 2 == 0:
        raise ValueError("Reference projection requires an odd smaller path")
    bead = np.arange(source)
    sample = np.arange(target)
    result = np.ones((target, source)) / source
    for frequency in range(1, (small + 1) // 2):
        difference = sample[:, None] / target - bead[None, :] / source
        result += 2 / source * np.cos(2 * np.pi * frequency * difference)
    return result


def polynomial(q, quadratic, quartic):
    """U=sum(a*x^2+b*x^4), using polynomial differentiation for F,H."""
    q = np.asarray(q)
    energy = np.zeros(q.shape[:-1])
    force, curvature = np.empty_like(q), np.empty_like(q)
    for coordinate, (a, b) in enumerate(zip(quadratic, quartic)):
        function = Polynomial([0, 0, a, 0, b])
        x = q[..., coordinate]
        energy += function(x)
        force[..., coordinate] = -function.deriv()(x)
        curvature[..., coordinate] = function.deriv(2)(x)
    return energy, force, curvature


def virials(q, force):
    """Coordinate virial only; complete when U has no explicit cell dependence."""
    return np.asarray(
        [np.triu(x.reshape(-1, 3).T @ f.reshape(-1, 3)) for x, f in zip(q, force)]
    )


def cell_strain_virial(energy, q, cell, step=1e-4):
    """-dE/dstrain at fixed fractional coordinates, in i-PI's convention.

    For component ij, q'=q@(I+eps*Eij) and h'=(I+eps*Eij).T@h.
    Both positions and cell change. `energy` must include explicit cell
    dependence and return either a scalar or one energy per bead.
    A five-point derivative has O(step^4) truncation error.
    """
    positions = np.asarray(q).reshape(len(q), -1, 3)
    value = np.asarray(energy(q, cell))
    result = np.zeros(value.shape + (3, 3))
    for i in range(3):
        for j in range(i, 3):
            direction = np.zeros((3, 3))
            direction[i, j] = 1
            values = []
            for displacement in [-2 * step, -step, step, 2 * step]:
                transform = np.eye(3) + displacement * direction
                moved = (positions @ transform).reshape(q.shape)
                values.append(energy(moved, transform.T @ cell))
            result[..., i, j] = -(
                values[0] - 8 * values[1] + 8 * values[2] - values[3]
            ) / (12 * step)
    return result


def volume_polynomial(q, cell, quadratic, quartic, coupling=0.0, volume0=120.0):
    """U=sum(a*(1+coupling*(V/V0-1))*x^2+b*x^4), including cell derivatives."""
    ratio = np.linalg.det(cell) / volume0
    a = np.asarray(quadratic)
    energy, force, curvature = polynomial(q, a * (1 + coupling * (ratio - 1)), quartic)
    # -V*dU/dV at fixed Cartesian positions supplements sum(q outer F).
    explicit = coupling * ratio * np.sum(a * q**2, axis=1)
    virial = virials(q, force) - explicit[:, None, None] * np.eye(3)
    return energy, force, curvature, virial


def sc_weights(count, alpha, temperature):
    """Published w_b,d_b; coefficient of sum(F^2/m) is w*d/omega_P^2."""
    weights = np.tile([2 / 3, 4 / 3], count // 2)
    double_commutator = np.tile([alpha / 6, (1 - alpha) / 12], count // 2)
    return weights, weights * double_commutator / (count * temperature) ** 2


def sc_reference(q, masses, quadratic, quartic, alpha, temperature):
    energy, force, curvature = polynomial(q, quadratic, quartic)
    base_virial = virials(q, force)
    count = len(q)
    w, c = sc_weights(count, alpha, temperature)
    mass = np.repeat(masses, 3)
    square = np.sum(force**2 / mass, axis=1)
    # -gradient(sum(F^2/m)) = 2*H*F/m for a diagonal polynomial Hessian.
    derivative = 2 * curvature * force / mass
    # -d(sum(F^2/m))/dstrain, not just -2*D(ordinary virial)[F/m].
    # With no explicit cell dependence, this is sum(q outer (-grad(F^2/m))).
    # Kapil et al., JCTC 15,3237 (2019), Eq. 9 / SI S1.26--S1.31.
    derivative_virial = virials(q, derivative)
    return dict(
        energy=energy,
        force=force,
        virial=base_virial,
        square=square,
        linear=w - 1,
        quadratic=c,
        derivative=derivative,
        derivative_virial=derivative_virial,
        correction=(w - 1) * energy + c * square,
        force1=(w - 1)[:, None] * force,
        force2=c[:, None] * derivative,
        virial1=(w - 1)[:, None, None] * base_virial,
        virial2=c[:, None, None] * derivative_virial,
    )


def spring_action(q, masses, temperature):
    difference = q - np.roll(q, 1, axis=0)
    return (
        0.5 * (len(q) * temperature) ** 2 * np.sum(difference**2 * np.repeat(masses, 3))
    )


def mode_mass_factors(count, temperature, style, frequencies):
    omega = 2 * count * temperature * np.sin(np.pi * np.arange(count) / count)
    factors = np.ones(count)
    if style == "manual":
        factors[1:] = (omega[1:] / frequencies) ** 2
    elif style == "pa-cmd":
        factors[1:] = (omega[1:] / frequencies[0]) ** 2
    elif style == "wmax-cmd":
        maximum, target = frequencies
        common = target * np.sqrt(1 + (maximum / omega[1]) ** 2)
        factors[1:] = (omega[1:] ** 2 + maximum**2) / common**2
    elif style != "rpmd":
        raise ValueError(style)
    return factors


def scaled_derivatives(q, quadratic, quartic, temperature):
    """Differentiate beta'*V(q_c+sqrt(beta'/beta)*dq; beta') analytically.

    Write s=beta'/beta. At s=1, dq/ds=dq/2 and d2q/ds2=-dq/4.
    """
    coordinates = q.shape[1]
    delta = q - q.mean(axis=0)
    energy, force, hessian = polynomial(q, quadratic, quartic)
    gradient = -force
    first = np.sum(gradient * delta / 2, axis=1)
    second = np.sum((hessian * delta**2 - gradient * delta) / 4, axis=1)
    e = np.mean(energy + first) + coordinates * temperature / 2
    ep = temperature * np.mean(2 * first + second)
    ep -= coordinates * temperature**2 / 2
    return np.array([e, ep])


def sc_scaled_derivatives(q, masses, quadratic, quartic, alpha, temperature):
    """Yamamoto Eqs.33-34 including explicit beta-dependence of SC action."""
    count, dimensions = q.shape
    delta = q - q.mean(axis=0)
    energy, force, curvature = polynomial(q, quadratic, quartic)
    gradient = -force
    mass = np.repeat(masses, 3)
    # Derivatives with respect to s=beta'/beta at s=1.
    g1 = curvature * delta / 2
    g2 = 24 * np.asarray(quartic) * q * delta**2 / 4 - curvature * delta / 4
    a = np.sum(gradient**2 / mass, axis=1)
    a1 = np.sum(2 * gradient * g1 / mass, axis=1)
    a2 = np.sum(2 * (g1**2 + gradient * g2) / mass, axis=1)
    u1 = np.sum(gradient * delta / 2, axis=1)
    u2 = np.sum((curvature * delta**2 - gradient * delta) / 4, axis=1)
    w, c = sc_weights(count, alpha, temperature)
    # d(s^3*A(s))/ds = 3A+A'; d2 = 6A+6A'+A''.
    e = dimensions * temperature / 2 + np.mean(w * (energy + u1) + c * (3 * a + a1))
    ep = -dimensions * temperature**2 / 2
    ep += temperature * np.mean(w * (2 * u1 + u2) + c * (6 * a + 6 * a1 + a2))
    return np.array([e, ep])


def permutation_reference(q, mass, temperature):
    """Explicit two-particle permutation sum, no recursion/exchange helpers."""
    count, atoms, _ = q.shape
    if atoms != 2:
        raise ValueError("The reference sums both permutations of two particles")
    energies, signs, longest = [], [], []
    common = np.sum(np.diff(q, axis=0) ** 2)
    for perm in itertools.permutations(range(atoms)):
        closure = np.sum((q[-1] - q[0, list(perm)]) ** 2)
        energies.append(0.5 * mass * (count * temperature) ** 2 * (common + closure))
        signs.append(1 if perm == (0, 1) else -1)
        longest.append(perm == (1, 0))
    weight = np.exp(-np.asarray(energies) / (count * temperature))
    probability = weight / weight.sum()
    return dict(
        distinct=probability[0],
        longest=atoms * np.dot(probability, longest),
        sign=np.dot(probability, signs),
        kinetic=1.5 * atoms * count * temperature
        - np.dot(probability, energies) / count,
    )
