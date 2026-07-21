"""Algorithms that deal with normal mode transformations and ring polymer contraction."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
from ipi.utils.messages import verbosity, info, warning
from ipi.utils.array_backend import xp, xp_setitem, backend

__all__ = [
    "nm_noop",
    "nm_trans",
    "nm_rescale",
    "nm_fft",
    "mk_nm_matrix",
    "mk_o_nm_matrix",
    "nm_eva",
    "o_nm_eva",
    "eco_eva",
    "eco_o_eva",
]


def mk_nm_matrix(nbeads):
    """
    Makes a matrix that transforms between the bead and normal mode
    representations.

    If we return from this function a matrix C, then we transform between the
    bead and normal mode representation using q_nm = C . q_b, q_b = C.T . q_nm

    Args:
       nbeads: The number of beads.
    """

    b2nm = np.zeros((nbeads, nbeads))
    b2nm[0, :] = np.sqrt(1.0)
    for j in range(nbeads):
        for i in range(1, nbeads // 2 + 1):
            b2nm[i, j] = np.sqrt(2.0) * np.cos(2 * np.pi * j * i / float(nbeads))
        for i in range(nbeads // 2 + 1, nbeads):
            b2nm[i, j] = np.sqrt(2.0) * np.sin(2 * np.pi * j * i / float(nbeads))
    if (nbeads % 2) == 0:
        b2nm[nbeads // 2, 0:nbeads:2] = 1.0
        b2nm[nbeads // 2, 1:nbeads:2] = -1.0
    return b2nm / np.sqrt(nbeads)


def nm_eva(nbeads):
    return 2 * np.array([np.sin(k * np.pi / nbeads) for k in range(nbeads)])


def o_nm_eva(nbeads):
    return 2 * np.array([np.sin(k * np.pi / (2 * nbeads)) for k in range(nbeads)])


def _eco_f(x):
    """Kernel f(x) = x^2 / ((x/2) coth(x/2) - 1) of the Eco objective function.

    Uses a Taylor expansion for small x to avoid numerical cancellation.
    """

    x = np.asarray(x, float)
    z = 0.5 * x
    small = z < 0.25
    z2 = np.where(small, z, 0.0) ** 2
    f_series = 4.0 / (1.0 / 3.0 - z2 / 45.0 + 2.0 * z2**2 / 945.0 - z2**3 / 4725.0)
    zb = np.where(small, 1.0, z)
    f_direct = np.where(small, 1.0, x) ** 2 / (zb / np.tanh(zb) - 1.0)
    return np.where(small, f_series, f_direct)


def _eco_fit(nbeads, xmax, y0=None):
    """Fits the dimensionless internal-mode parameters y_k = beta*hbar*omega_k
    of the Eco path integral, minimizing the rms fractional error in the
    radius of gyration of harmonic oscillators with 0 <= beta*hbar*omega <= xmax.
    Follows the reference implementation in the supplementary material of
    Zeng & Manolopoulos, "Economised path integrals", arXiv:2607.06414
    (https://arxiv.org/abs/2607.06414): safe Newton iterations
    with an eigenvalue-shifted Hessian and a line search that keeps the y_k
    positive and in ascending order. Starts from the Matsubara frequencies,
    or from the initial guess y0 (e.g. a previous solution) if given.

    Returns an array of nbeads//2 optimized y_k, in ascending order.
    """

    nfree = nbeads // 2
    # each internal mode is doubly degenerate (y_k = y_{nbeads-k}) except
    # for the middle one when nbeads is even
    mult = np.full(nfree, 2.0)
    if nbeads % 2 == 0:
        mult[-1] = 1.0

    # midpoint grid of oscillator frequencies in (0, xmax), 10 per unit x
    m = max(round(10 * xmax), 100)
    x = (np.arange(m) + 0.5) * (xmax / m)
    f = _eco_f(x)
    x2 = x[:, np.newaxis] ** 2

    def objfun(y):
        """s(y) = (1/2m) sum_j r_j^2, with its gradient and Hessian."""
        d = 1.0 / (y**2 + x2)
        e = f[:, np.newaxis] * mult * d
        r = e.sum(axis=1) - 1.0
        dg = -2.0 * d * e * y
        d2 = 2.0 * d * d * e * (3.0 * y**2 - x2)
        s = 0.5 * (r**2).mean()
        g = (r[:, np.newaxis] * dg).mean(axis=0)
        h = (dg.T @ dg + np.diag(r @ d2)) / m
        return s, g, h

    if y0 is not None:
        y = np.array(y0, float)
    else:
        y = 2.0 * np.pi * np.arange(1, nfree + 1, dtype=float)  # Matsubara guess
    s, g, h = objfun(y)
    for _ in range(500):
        # Newton shift, offsetting the Hessian eigenvalues to get a descent direction
        eva, vec = np.linalg.eigh(h)
        delta = max(1e-16 * eva[-1], -2.0 * eva[0])
        dy = vec @ (-(vec.T @ g) / (eva + delta))
        # backtracking line search that preserves ordering and positivity
        sp = s
        c = 1.0
        for _ in range(60):
            z = y + c * dy
            if z[0] >= 0.0 and np.all(np.diff(z) >= 0.0):
                sn, gn, hn = objfun(z)
                if sn <= s:
                    y, s, g, h = z, sn, gn, hn
                    break
            c *= 0.5
        else:
            break
        # stops when the relative decrease of the objective becomes negligible
        if sp - s <= 1e-12 * sp:
            break
    else:
        # when many modes fit an easy target the minimum is a flat valley and
        # the loop can spend all iterations shaving negligible amounts off an
        # already-excellent fit; an error is raised only if the exhausted
        # optimisation is still far from a stationary point (large gradient)
        # of a good fit (rms error in R^2 above ~1e-4, i.e. sqrt(2e-8))
        if s > 1e-8 and np.abs(g).max() > 1e-6:
            raise ValueError(
                "Eco frequency optimisation did not converge in 500 iterations for "
                "nbeads=%d, xmax=%g (rms fractional error in R^2 = %g); check that the "
                "maximum frequency and the temperature are physically sensible."
                % (nbeads, xmax, np.sqrt(2.0 * s))
            )

    info(
        " @nmtransform: Eco fit for nbeads=%d, xmax=%g: rms fractional error in R^2 = %g"
        % (nbeads, xmax, np.sqrt(2.0 * s)),
        verbosity.medium,
    )
    return y


def eco_eva(nbeads, xmax, y0=None):
    """Computes dimensionless eigenvalues of the Eco ring-polymer springs,
    optimized to reproduce the radii of gyration of harmonic oscillators
    with frequencies 0 <= beta*hbar*omega <= xmax. Defined so that
    omega_k = omegan * eco_eva(nbeads, xmax)_k, in analogy with nm_eva.
    An initial guess y0 for the nbeads//2 free parameters (e.g. the solution
    at a nearby temperature) can be given to speed up the fit.
    See Zeng & Manolopoulos, "Economised path integrals", arXiv:2607.06414
    (https://arxiv.org/abs/2607.06414).
    """

    if xmax <= 0:
        raise ValueError("Eco path integrals require a positive maximum frequency.")
    if nbeads == 1:
        return np.zeros(1)
    # the guess must satisfy the line-search invariants, else start from scratch
    if y0 is not None and not (
        len(y0) == nbeads // 2 and np.all(y0 > 0) and np.all(np.diff(y0) >= 0)
    ):
        y0 = None
    y = _eco_fit(nbeads, float(xmax), y0)
    eva = np.zeros(nbeads)
    for k in range(1, nbeads):
        eva[k] = y[min(k, nbeads - k) - 1]
    return eva / nbeads


def _eco_open_t(x):
    """End-to-end target t(x) = (2/x) tanh(x/2): the exact variance of the
    quantum end-to-end distribution N(Delta) of a harmonic oscillator with
    x = beta*hbar*omega, in units of beta*hbar^2/m.

    Uses a Taylor expansion for small x to avoid numerical cancellation.
    """

    x = np.asarray(x, float)
    small = x < 0.5
    z2 = np.where(small, 0.5 * x, 0.0) ** 2
    t_series = 1.0 - z2 / 3.0 + 2.0 * z2**2 / 15.0 - 17.0 * z2**3 / 315.0
    xb = np.where(small, 1.0, x)
    return np.where(small, t_series, 2.0 / xb * np.tanh(xb / 2.0))


def _eco_open_g(x):
    """Open-path gyration target g(x) = (x coth x - 1)/(2 x^2): the exact
    radius of gyration of a free-ended (open) harmonic path with
    x = beta*hbar*omega, in units of beta*hbar^2/m.

    Uses a Taylor expansion for small x to avoid numerical cancellation.
    """

    x = np.asarray(x, float)
    x2 = np.where(x < 0.5, x, 0.0) ** 2
    g_series = 1.0 / 6.0 - x2 / 90.0 + x2**2 / 945.0 - x2**3 / 9450.0
    xb = np.where(x < 0.5, 1.0, x)
    return np.where(x < 0.5, g_series, (xb / np.tanh(xb) - 1.0) / (2.0 * xb**2))


def _eco_open_objfun(wt, x2, r0):
    """Builds a least-squares objective for residuals of the generic form
    r_j = sum_k wt_jk / (y_k^2 + x_j^2) + r0_j, returning value, gradient
    and Hessian as in the closed-path _eco_fit."""

    def objfun(y):
        d = 1.0 / (y**2 + x2)
        e = wt * d
        r = e.sum(axis=1) + r0
        dg = -2.0 * d * e * y
        d2 = 2.0 * d * d * e * (3.0 * y**2 - x2)
        s = 0.5 * (r**2).mean()
        g = (r[:, np.newaxis] * dg).mean(axis=0)
        h = (dg.T @ dg + np.diag(r @ d2)) / len(r)
        return s, g, h

    return objfun


def _eco_open_newton(y, objfun, groups, ycap, check=False, label=""):
    """Safe Newton minimization with an eigenvalue-shifted Hessian and a
    backtracking line search that keeps the y_k positive, below ycap, and,
    within each index group in `groups`, in ascending order (odd- and even-k
    open-path modes play different roles and need not be mutually ordered).
    The cap prevents redundant modes, whose contribution to the objective is
    negligible, from drifting towards infinite spring stiffness."""

    s, g, h = objfun(y)
    for _ in range(500):
        eva, vec = np.linalg.eigh(h)
        delta = max(1e-16 * eva[-1], -2.0 * eva[0])
        dy = vec @ (-(vec.T @ g) / (eva + delta))
        sp = s
        c = 1.0
        for _ in range(60):
            z = y + c * dy
            if (
                z.min() >= 0.0
                and z.max() <= ycap
                and all(np.all(np.diff(z[gr]) >= 0.0) for gr in groups)
            ):
                sn, gn, hn = objfun(z)
                if sn <= s:
                    y, s, g, h = z, sn, gn, hn
                    break
            c *= 0.5
        else:
            break
        # stops when the relative decrease of the objective becomes negligible
        if sp - s <= 1e-12 * sp:
            break
    else:
        # see the analogous guard in _eco_fit: an error is raised only if the
        # exhausted optimisation is far from a stationary point of a good fit
        # (thresholds looser than the closed-path ones, as the open objective
        # retains an irreducible residual when nbeads is small)
        if check and s > 1e-4 and np.abs(g).max() > 1e-5:
            raise ValueError(
                "Eco open-path frequency optimisation did not converge in 500 "
                "iterations%s (rms residual = %g); check that the maximum "
                "frequency and the temperature are physically sensible."
                % (label, np.sqrt(2.0 * s))
            )
    return y, s


def _eco_open_fit(nbeads, xmax, y0=None):
    """Fits the dimensionless mode parameters y_k = beta*hbar*omega_k,
    k = 1..nbeads-1, of the Eco open path integral. Minimizes jointly the
    rms fractional errors, over harmonic oscillators with
    0 <= x = beta*hbar*omega <= xmax, of
      (A) the end-to-end variance including the endpoint-kernel term 1/P:
          8 sum_{k odd} cos^2(k pi/2P)/(y_k^2+x^2) + 1/P  vs  (2/x)tanh(x/2)
      (B) the open-path radius of gyration:
          sum_k 1/(y_k^2+x^2)  vs  (x coth x - 1)/(2x^2).
    Only odd-k modes couple to the end-to-end vector, so the joint optimum is
    found robustly by first fitting the odd modes on (A) and the even modes on
    the remainder of (B), then polishing with Newton steps on the joint
    objective. A warm-start y0 (e.g. the solution at a nearby temperature)
    skips the two construction stages.

    Returns an array of nbeads-1 optimized y_k (ascending within each parity
    class, but not necessarily overall).
    """

    P = nbeads
    m = max(round(10 * xmax), 100)
    x = (np.arange(m) + 0.5) * (xmax / m)
    t = _eco_open_t(x)
    gy = _eco_open_g(x)
    x2 = x[:, np.newaxis] ** 2

    # the end-to-end variance of oscillators with x >~ P falls below the
    # endpoint-kernel variance 1/P and cannot be represented by the truncated
    # chain, so block A is fitted only where the kernel is at most half the
    # exact variance; chasing the unrepresentable tail would push the odd
    # frequencies to infinity
    in_a = t >= 2.0 / P
    if not in_a.all():
        warning(
            "Eco open path: nbeads=%d is too small to represent the end-to-end "
            "distribution up to beta*hbar*omega_max=%g; its fit is restricted "
            "to beta*hbar*omega < %g."
            % (nbeads, xmax, x[in_a].max() if in_a.any() else 0.0),
            verbosity.low,
        )

    k = np.arange(1, P)
    odd = k % 2 == 1
    wa = np.where(odd, 8.0 * np.cos(k * np.pi / (2 * P)) ** 2, 0.0)
    # block A: r = [sum_k wa_k d_k + 1/P]/t - 1;  block B: r = [sum_k d_k]/g - 1
    wt_a = np.outer(1.0 / t[in_a], wa)
    x2_a = x2[in_a]
    r0_a = 1.0 / (P * t[in_a]) - 1.0
    wt_b = np.outer(1.0 / gy, np.ones(P - 1))
    groups = [np.where(odd)[0], np.where(~odd)[0]]

    # frequencies above twice the stiffest Trotter mode (y ~ 2P) contribute
    # negligibly to either observable and are only detrimental to integration
    ycap = 4.0 * P

    if y0 is not None:
        y = np.array(y0, float)
    else:
        y = np.pi * k.astype(float)  # open-path Matsubara guess
        # stage 1: odd modes alone on the end-to-end block
        if in_a.any():
            yo, _ = _eco_open_newton(
                y[odd],
                _eco_open_objfun(wt_a[:, odd], x2_a, r0_a),
                [np.arange(odd.sum())],
                ycap,
            )
            y[odd] = yo
        # stage 2: even modes alone on the remainder of the gyration block
        rem = gy - (1.0 / (y[odd][np.newaxis, :] ** 2 + x2)).sum(axis=1)
        rem = np.maximum(rem, 1e-6 * gy)
        if (~odd).any():
            ye, _ = _eco_open_newton(
                y[~odd],
                _eco_open_objfun(
                    np.outer(1.0 / rem, np.ones((~odd).sum())), x2, -np.ones(m)
                ),
                [np.arange((~odd).sum())],
                ycap,
            )
            y[~odd] = ye

    # joint polish of both blocks
    y, s = _eco_open_newton(
        y,
        _eco_open_objfun(
            np.vstack([wt_a, wt_b]),
            np.concatenate([x2_a, x2]),
            np.concatenate([r0_a, -np.ones(m)]),
        ),
        groups,
        ycap,
        check=True,
        label=" for nbeads=%d, xmax=%g" % (nbeads, xmax),
    )

    # separate rms fractional errors of the two fitted observables (block A
    # reported over its representable window)
    d = 1.0 / (y[np.newaxis, :] ** 2 + x2)
    rms_a = (
        np.sqrt((((((d * wa).sum(axis=1) + 1.0 / P) / t - 1.0) ** 2)[in_a]).mean())
        if in_a.any()
        else float("nan")
    )
    rms_b = np.sqrt(((d.sum(axis=1) / gy - 1.0) ** 2).mean())
    info(
        " @nmtransform: Eco open-path fit for nbeads=%d, xmax=%g: rms fractional "
        "error in end-to-end variance = %g, in gyration = %g"
        % (nbeads, xmax, rms_a, rms_b),
        verbosity.medium,
    )
    return y


def eco_o_eva(nbeads, xmax, y0=None):
    """Computes dimensionless eigenvalues of the Eco springs of an open path,
    optimized to reproduce the end-to-end distribution (hence the particle
    momentum distribution) and the radius of gyration of open harmonic paths
    with frequencies 0 <= beta*hbar*omega <= xmax. Defined so that
    omega_k = omegan * eco_o_eva(nbeads, xmax)_k, in analogy with o_nm_eva.
    An initial guess y0 for the nbeads-1 free parameters (e.g. the solution
    at a nearby temperature) can be given to speed up the fit.
    """

    if xmax <= 0:
        raise ValueError("Eco path integrals require a positive maximum frequency.")
    if nbeads == 1:
        return np.zeros(1)
    # the guess must satisfy the line-search invariants (positivity, cap,
    # ordering within each parity class), else start from scratch
    if y0 is not None:
        y0 = np.asarray(y0, float)
        if not (
            len(y0) == nbeads - 1
            and np.all(y0 > 0)
            and np.all(y0 <= 4.0 * nbeads)
            and np.all(np.diff(y0[0::2]) >= 0)
            and np.all(np.diff(y0[1::2]) >= 0)
        ):
            y0 = None
    y = _eco_open_fit(nbeads, float(xmax), y0)
    eva = np.zeros(nbeads)
    eva[1:] = y
    return eva / nbeads


def mk_o_nm_matrix(nbeads):
    """
    Makes a matrix that transforms between the bead and the (open path) normal mode
    representations.
    """
    # here define the orthogonal transformation matrix for the open path
    b2o_nm = np.zeros((nbeads, nbeads))
    b2o_nm[0, :] = np.sqrt(1.0)
    for j in range(0, nbeads):
        for i in range(1, nbeads):
            b2o_nm[i, j] = np.sqrt(2.0) * np.cos(np.pi * (j + 0.5) * i / float(nbeads))
    return b2o_nm / np.sqrt(nbeads)


def mk_rs_matrix(nb1, nb2):
    """Makes a matrix that transforms a path with `nb1` beads to one with `nb2` beads.

    If we return from this function a matrix T, then we transform between the
    system with nb1 bead and the system of nb2 beads using q_2 = T . q_1

    Args:
       nb1: The initial number of beads.
       nb2: The final number of beads.
    """

    if nb1 == nb2:
        return np.identity(nb1, float)
    elif nb1 > nb2:
        b1_nm = mk_nm_matrix(nb1)
        nm_b2 = mk_nm_matrix(nb2).T

        # builds the "reduction" matrix that picks the normal modes we want to keep
        b1_b2 = np.zeros((nb2, nb1), float)
        b1_b2[0, 0] = 1.0
        for i in range(1, nb2 // 2 + 1):
            b1_b2[i, i] = 1.0
            b1_b2[nb2 - i, nb1 - i] = 1.0
        if nb2 % 2 == 0:
            # if we are contracting down to an even number of beads, then we have to
            # pick just one of the last degenerate modes to match onto the single
            # stiffest mode in the new path
            b1_b2[nb2 // 2, nb1 - nb2 // 2] = 0.0

        rs_b1_b2 = np.dot(nm_b2, np.dot(b1_b2, b1_nm))
        return rs_b1_b2 * np.sqrt(float(nb2) / float(nb1))
    else:
        return mk_rs_matrix(nb2, nb1).T * (float(nb2) / float(nb1))


def mk_o_rs_matrix(nb1, nb2):
    """Makes a matrix that transforms a path with `nb1` beads to one with `nb2` beads.

    If we return from this function a matrix T, then we transform between the
    system with nb1 bead and the system of nb2 beads using q_2 = T . q_1

    Args:
       nb1: The initial number of beads.
       nb2: The final number of beads.
    """

    if nb1 == nb2:
        return np.identity(nb1, float)
    elif nb1 > nb2:
        b1_nm = mk_o_nm_matrix(nb1)
        nm_b2 = mk_o_nm_matrix(nb2).T

        # builds the "reduction" matrix that picks the normal modes we want to keep
        b1_b2 = np.zeros((nb2, nb1), float)
        b1_b2[0, 0] = 1.0
        for i in range(1, nb2 // 2 + 1):
            b1_b2[i, i] = 1.0
            b1_b2[nb2 - i, nb1 - i] = 1.0
        if nb2 % 2 == 0:
            # if we are contracting down to an even number of beads, then we have to
            # pick just one of the last degenerate modes to match onto the single
            # stiffest mode in the new path
            b1_b2[nb2 // 2, nb1 - nb2 // 2] = 0.0

        rs_b1_b2 = np.dot(nm_b2, np.dot(b1_b2, b1_nm))
        return rs_b1_b2 * np.sqrt(float(nb2) / float(nb1))
    else:
        return mk_o_rs_matrix(nb2, nb1).T * (float(nb2) / float(nb1))


class nm_noop(object):
    """A no-op NM transformation for classical trajectories"""

    def __init__(self, nbeads, open_paths=None):
        """Initializes nm_noop."""
        if nbeads > 1:
            raise ValueError(
                "Shouldn't use a noop transformation for ring-polymer systems "
            )

    def b2nm(self, q):
        return q

    def nm2b(self, qnm):
        return qnm


def _apply(matrix, q):
    """Multiplies `matrix` by `q`, contracting the first axis of `q`, that can
    have any number of dimensions."""

    if q.ndim <= 2:
        return matrix @ q
    q2 = q.reshape(q.shape[0], -1)
    return (matrix @ q2).reshape(matrix.shape[0], *q.shape[1:])


def open_columns(open_paths):
    """Indices of the coordinates of the atoms with open paths, as an array of
    the backend."""

    cols = np.array([3 * io + k for io in open_paths for k in range(3)], int)
    return xp.asarray(cols)


class nm_trans(object):
    """Uses matrix multiplication to do normal mode transformations.

    Attributes:
       _b2nm: The matrix to transform between the bead and normal mode
          representations.
       _nm2b: The matrix to transform between the normal mode and bead
          representations.
    """

    def __init__(self, nbeads, open_paths=None):
        """Initializes nm_trans.

        Args:
           nbeads: The number of beads.
        """

        self._b2nm = xp.asarray(mk_nm_matrix(nbeads))
        self._nm2b = self._b2nm.T
        if open_paths is None:
            open_paths = []
        self._open = open_paths
        self._open_cols = open_columns(open_paths)
        # definition of the transformation also with the open path matrx
        self._b2o_nm = xp.asarray(mk_o_nm_matrix(nbeads))
        self._o_nm2b = self._b2o_nm.T

    def b2nm(self, q):
        """Transforms a matrix to the normal mode representation.

        Args:
           q: A matrix with nbeads rows, in the bead representation.
        """

        qnm = self._b2nm @ q
        if len(self._open) > 0:
            # the atoms with open paths use the open-path transformation
            cols = self._open_cols
            qnm = xp_setitem(qnm, (slice(None), cols), self._b2o_nm @ q[:, cols])
        return qnm

    def nm2b(self, qnm):
        """Transforms a matrix to the bead representation.

        Args:
           q: A matrix with nbeads rows, in the normal mode representation.
        """

        q = self._nm2b @ qnm
        if len(self._open) > 0:
            cols = self._open_cols
            q = xp_setitem(q, (slice(None), cols), self._o_nm2b @ qnm[:, cols])
        return q


class nm_rescale(object):
    """Uses matrix multiplication to do ring polymer contraction or expansion
    between different numbers of beads.

    Attributes:
       _b1tob2: The matrix to transform between a ring polymer with 'nbeads1'
          beads and another with 'nbeads2' beads.
       _b2tob1: The matrix to transform between a ring polymer with 'nbeads2'
          beads and another with 'nbeads1' beads.
    """

    def __init__(self, nbeads1, nbeads2, open_paths=None):
        """Initializes nm_rescale.

        Args:
           nbeads1: The initial number of beads.
           nbeads2: The rescaled number of beads.
        """

        self.nbeads1 = nbeads1
        self.nbeads2 = nbeads2

        if open_paths is None:
            open_paths = []
        self._open = open_paths
        self._open_cols = open_columns(open_paths)

        self.noop = nbeads1 == nbeads2
        if self.noop:
            # a contraction between the same numbers of beads only copies
            self._b1tob2 = self._b2tob1 = self._o_b1tob2 = self._o_b2tob1 = None
        else:
            b1tob2 = mk_rs_matrix(nbeads1, nbeads2)
            self._b1tob2 = xp.asarray(b1tob2)
            self._b2tob1 = xp.asarray(b1tob2.T * (float(nbeads1) / float(nbeads2)))
            # definition of the scaling also using the open case normal mode matrix transformations
            o_b1tob2 = mk_o_rs_matrix(nbeads1, nbeads2)
            self._o_b1tob2 = xp.asarray(o_b1tob2)
            self._o_b2tob1 = xp.asarray(o_b1tob2.T * (float(nbeads1) / float(nbeads2)))

    def _rescale(self, q, matrix, o_matrix):
        """Applies the rescaling `matrix` to `q`, and `o_matrix` to the atoms
        with open paths."""

        # inputs can also be numpy arrays, e.g. the extras returned by forcefields
        q = xp.asarray(q)
        if self.noop:
            # still must return a copy, as the contraction is meant to return new data, not a view
            return xp.asarray(q, copy=True)
        if len(self._open) == 0:
            return _apply(matrix, q)
        if q.ndim != 2:
            # this applies the open path contraction to EVERYTHING because otherwise we don't know how to handle
            # the fact that only some beads are open. clearly this is a hack, and in practice the point is that
            # a "per bead" NM transformation of the potential is not well-defined when different beads have different
            # NM transformations
            return _apply(o_matrix, q)
        # this applies to bead vector properties (e.g. positions, forces)
        q_scal = matrix @ q
        cols = self._open_cols
        return xp_setitem(q_scal, (slice(None), cols), o_matrix @ q[:, cols])

    def b1tob2(self, q):
        """Transforms a matrix from one value of beads to another.

        Args:
           q: A matrix with nbeads1 rows, in the bead representation.
        """

        return self._rescale(q, self._b1tob2, self._o_b1tob2)

    def b2tob1(self, q):
        """Transforms a matrix from one value of beads to another.

        Args:
           q: A matrix with nbeads2 rows, in the bead representation.
        """

        return self._rescale(q, self._b2tob1, self._o_b2tob1)


class nm_fft(object):
    """Uses Fast Fourier transforms to do normal mode transformations.

    The real FFT of the bead coordinates along the beads axis gives the
    cosine (real part) and sine (imaginary part) components of the normal
    modes, which are scaled to make the transformation orthogonal and arranged
    in the same order as the rows of `mk_nm_matrix`. The atoms with open paths
    are transformed with the open-path matrix. With the numpy backend, PyFFTW
    is used if it is available.

    Attributes:
       nbeads: The number of beads.
       natoms: The number of atoms.
       nmodes: The number of cosine components beyond the centroid.
       nsin: The number of sine components.
    """

    def __init__(
        self, nbeads, natoms, open_paths=None, n_threads=1, single_precision=False
    ):
        """Initializes nm_fft.

        Args:
           nbeads: The number of beads.
           natoms: The number of atoms.
           open_paths: Indices of the atoms with open paths.
           n_threads: Number of threads used by PyFFTW.
           single_precision: Whether the FFT is computed in single precision.
        """

        self.nbeads = nbeads
        self.natoms = natoms
        self.n_threads = n_threads
        self.single_precision = single_precision
        self.nmodes = nbeads // 2
        self.nsin = nbeads - self.nmodes - 1
        # scaling of the FFT components that makes the transformation orthogonal
        weights = np.full(self.nmodes + 1, np.sqrt(2.0))
        weights[0] = 1.0
        if nbeads % 2 == 0:
            weights[-1] = 1.0
        self._weights = xp.asarray(weights[:, np.newaxis] / np.sqrt(nbeads))
        if open_paths is None:
            open_paths = []
        self._open = open_paths
        self._open_cols = open_columns(open_paths)
        # for atoms with open path we still use the matrix transformation
        self._b2o_nm = xp.asarray(mk_o_nm_matrix(nbeads))
        self._o_nm2b = self._b2o_nm.T
        self._setup_fft()

    def _setup_fft(self):
        """Defines `_rfft(q)`, the real FFT along the beads axis, and
        `_irfft(re, im)`, its inverse given the real parts of the FFT
        components and the imaginary parts of the components 1 to nsin
        (the others being zero)."""

        nbeads, natoms = self.nbeads, self.natoms
        if backend() == "numpy":
            try:
                import pyfftw

                info("Import of PyFFTW successful", verbosity.medium)
                qbuf = pyfftw.empty_aligned(
                    (nbeads, 3 * natoms),
                    dtype="float32" if self.single_precision else "float64",
                    n=16,
                )
                cbuf = pyfftw.empty_aligned(
                    (nbeads // 2 + 1, 3 * natoms),
                    dtype="complex64" if self.single_precision else "complex128",
                    n=16,
                )
                pyfftw.config.NUM_THREADS = self.n_threads
                forward = pyfftw.FFTW(
                    qbuf,
                    cbuf,
                    axes=(0,),
                    direction="FFTW_FORWARD",
                    threads=self.n_threads,
                )
                backward = pyfftw.FFTW(
                    cbuf,
                    qbuf,
                    axes=(0,),
                    direction="FFTW_BACKWARD",
                    threads=self.n_threads,
                )

                nsin = self.nsin

                def rfft(q):
                    qbuf[:] = q
                    forward()
                    return cbuf

                def irfft(re, im):
                    # writes directly in the real and imaginary parts of the buffer
                    cbuf.real[:] = re
                    cbuf.imag[0] = 0.0
                    cbuf.imag[1 : 1 + nsin] = im
                    cbuf.imag[1 + nsin :] = 0.0
                    backward()
                    return qbuf.copy()

                self._rfft, self._irfft = rfft, irfft
                return
            except ImportError:
                info(
                    "Import of PyFFTW unsuccessful, using NumPy library instead. ",
                    verbosity.low,
                )
                if nbeads > 1:
                    warning(
                        "Install PyFFTW to avoid slow normal-modes integration",
                        verbosity.low,
                    )

        nsin = self.nsin

        def to_complex(re, im):
            # the imaginary part is zero for the centroid and, for an even
            # number of beads, for the last component
            im = xp.concat(
                [xp.zeros_like(re[:1]), im, xp.zeros_like(re[1 + nsin :])], axis=0
            )
            return re + 1j * im

        if self.single_precision:

            def rfft(q):
                return xp.fft.rfft(xp.astype(q, xp.float32), axis=0)

            def irfft(re, im):
                c = xp.astype(to_complex(re, im), xp.complex64)
                return xp.astype(xp.fft.irfft(c, n=nbeads, axis=0), xp.float64)

        else:

            def rfft(q):
                return xp.fft.rfft(q, axis=0)

            def irfft(re, im):
                return xp.fft.irfft(to_complex(re, im), n=nbeads, axis=0)

        self._rfft, self._irfft = rfft, irfft

    def b2nm(self, q):
        """Transforms a matrix to the normal mode representation.

        Args:
           q: A matrix with nbeads rows and 3*natoms columns,
              in the bead representation.
        """

        if self.nbeads == 1:
            return q
        c = self._rfft(q) * self._weights
        # cosine components first, then sine components in reverse order
        qnm = xp.concat(
            [xp.real(c), xp.flip(xp.imag(c[1 : 1 + self.nsin]), axis=0)], axis=0
        )
        if self.single_precision:
            qnm = xp.astype(qnm, q.dtype)
        if len(self._open) > 0:
            cols = self._open_cols
            qnm = xp_setitem(qnm, (slice(None), cols), self._b2o_nm @ q[:, cols])
        return qnm

    def nm2b(self, qnm):
        """Transforms a matrix to the bead representation.

        Args:
           qnm: A matrix with nbeads rows and 3*natoms columns,
              in the normal mode representation.
        """

        if self.nbeads == 1:
            return qnm
        nmodes = self.nmodes
        # real (cosine) and imaginary (sine, stored in reverse order) parts
        # of the non-zero FFT components
        re = qnm[: nmodes + 1] / self._weights
        im = xp.flip(qnm[nmodes + 1 :], axis=0) / self._weights[1 : 1 + self.nsin]
        q = self._irfft(re, im)
        if len(self._open) > 0:
            cols = self._open_cols
            q = xp_setitem(q, (slice(None), cols), self._o_nm2b @ qnm[:, cols])
        return q
