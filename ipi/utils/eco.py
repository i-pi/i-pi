"""Economised (Eco) path integral spring frequencies, for closed and open
paths: least-squares fits of the free-ring-polymer normal-mode spectrum to
exact harmonic-oscillator observables, following Zeng & Manolopoulos,
"Economised path integrals", arXiv:2607.06414
(https://arxiv.org/abs/2607.06414), and its generalization to the open
(truncated-chain) paths of Kapil, Cuzzocrea & Ceriotti, JPCB 122, 6048
(2018) used for momentum distribution estimation."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.messages import verbosity, info, warning

__all__ = [
    "eco_eva",
    "eco_o_eva",
    "eco_o_kernel",
]


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
    k = 1..nbeads-1, of the Eco open path integral, using the analytical
    endpoint-kernel variance of eco_o_kernel. Minimizes jointly the rms
    fractional errors, over harmonic oscillators with
    0 <= x = beta*hbar*omega <= xmax, of
      (A) the end-to-end variance including the endpoint kernel:
          8 sum_{k odd} cos^2(k pi/2P)/(y_k^2+x^2) + u  vs  (2/x)tanh(x/2)
      (B) the open-path radius of gyration:
          sum_k 1/(y_k^2+x^2)  vs  (x coth x - 1)/(2x^2).
    Only odd-k modes couple to the end-to-end vector, so the joint optimum
    is found robustly by first fitting the odd modes on (A) and the even
    modes on the remainder of (B), then polishing with Newton steps on the
    joint objective. A warm-start y0 (e.g. the solution at a nearby
    temperature) skips the two construction stages.

    Returns (y, u): the nbeads-1 optimized y_k (ascending within each
    parity class, but not necessarily overall) and the kernel variance u.
    """

    P = nbeads
    u = eco_o_kernel(nbeads, xmax)
    m = max(round(10 * xmax), 100)
    x = (np.arange(m) + 0.5) * (xmax / m)
    t = _eco_open_t(x)
    gy = _eco_open_g(x)
    x2 = x[:, np.newaxis] ** 2

    k = np.arange(1, P)
    odd = k % 2 == 1
    wa = np.where(odd, 8.0 * np.cos(k * np.pi / (2 * P)) ** 2, 0.0)
    # block A: r = [sum_k wa_k d_k + u]/t - 1;  block B: r = [sum_k d_k]/g - 1
    wt_a = np.outer(1.0 / t, wa)
    r0_a = u / t - 1.0
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
        yo, _ = _eco_open_newton(
            y[odd],
            _eco_open_objfun(wt_a[:, odd], x2, r0_a),
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
            np.concatenate([x2, x2]),
            np.concatenate([r0_a, -np.ones(m)]),
        ),
        groups,
        ycap,
        check=True,
        label=" for nbeads=%d, xmax=%g" % (nbeads, xmax),
    )

    # separate rms fractional errors of the two fitted observables
    d = 1.0 / (y[np.newaxis, :] ** 2 + x2)
    rms_a = np.sqrt(((((d * wa).sum(axis=1) + u) / t - 1.0) ** 2).mean())
    rms_b = np.sqrt(((d.sum(axis=1) / gy - 1.0) ** 2).mean())
    info(
        " @nmtransform: Eco open-path fit for nbeads=%d, xmax=%g: rms fractional "
        "error in end-to-end variance = %g, in gyration = %g; kernel variance "
        "u = %g (Trotter: %g)" % (nbeads, xmax, rms_a, rms_b, u, 1.0 / P),
        verbosity.medium,
    )
    if u < 0.99 / P:
        warning(
            "Eco open path: the endpoint kernel was economised to u = %g "
            "(Trotter: 1/P = %g) because nbeads=%d < beta*hbar*omega_max = %g. "
            "End-to-end (momentum) estimators must smear with the economised "
            "kernel variance u*beta*hbar^2/m (e.g. get_np_rad.py -wmax)."
            % (u, 1.0 / P, nbeads, xmax),
            verbosity.low,
        )
    return y, u


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
    y, _ = _eco_open_fit(nbeads, float(xmax), y0)
    eva = np.zeros(nbeads)
    eva[1:] = y
    return eva / nbeads


def eco_o_kernel(nbeads, xmax):
    """Computes the economised endpoint-kernel variance of an Eco open path,
    in units of beta*hbar^2/m,

       u = (2/xmax) tanh(xmax / (2 nbeads)),

    the smearing that the two endpoint half-segments accumulate for a
    harmonic oscillator at the maximum fitted frequency (the free-particle
    result 1/nbeads of the Trotter factorization, recovered for
    nbeads >> xmax, would exceed the exact end-to-end variance of stiff
    oscillators when nbeads < xmax, making them unrepresentable).
    End-to-end (momentum distribution) estimators must smear the sampled
    end-to-end vectors with a Gaussian of this variance.
    """

    if xmax <= 0:
        raise ValueError("Eco path integrals require a positive maximum frequency.")
    if xmax < 1e-6:
        return 1.0 / nbeads
    return 2.0 * np.tanh(0.5 * xmax / nbeads) / xmax
