"""Algorithms that deal with normal mode transformations and ring polymer contraction."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.depend import dstrip
from ipi.utils.messages import verbosity, info, warning

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

        self._b2nm = mk_nm_matrix(nbeads)
        self._nm2b = self._b2nm.T
        if open_paths is None:
            open_paths = []
        self._open = open_paths
        # definition of the transformation also with the open path matrx
        self._b2o_nm = mk_o_nm_matrix(nbeads)
        self._o_nm2b = self._b2o_nm.T

    def b2nm(self, q):
        """Transforms a matrix to the normal mode representation.

        Args:
           q: A matrix with nbeads rows, in the bead representation.
        """

        qnm = np.tensordot(self._b2nm, q, axes=(1, 0))
        if len(self._open) > 0:
            for (
                io
            ) in (
                self._open
            ):  # does separately the transformation for the atom that are marked as open paths
                qnm[:, 3 * io] = np.tensordot(self._b2o_nm, q[:, 3 * io], axes=(1, 0))
                qnm[:, 3 * io + 1] = np.tensordot(
                    self._b2o_nm, q[:, 3 * io + 1], axes=(1, 0)
                )
                qnm[:, 3 * io + 2] = np.tensordot(
                    self._b2o_nm, q[:, 3 * io + 2], axes=(1, 0)
                )

        return qnm

    def nm2b(self, qnm):
        """Transforms a matrix to the bead representation.

        Args:
           q: A matrix with nbeads rows, in the normal mode representation.
        """

        q = np.tensordot(self._nm2b, qnm, axes=(1, 0))
        if len(self._open) > 0:
            for io in self._open:
                q[:, 3 * io] = np.tensordot(self._o_nm2b, qnm[:, 3 * io], axes=(1, 0))
                q[:, 3 * io + 1] = np.tensordot(
                    self._o_nm2b, qnm[:, 3 * io + 1], axes=(1, 0)
                )
                q[:, 3 * io + 2] = np.tensordot(
                    self._o_nm2b, qnm[:, 3 * io + 2], axes=(1, 0)
                )

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

        if nbeads1 == nbeads2:
            # avoids doing much if we try to contract between sambe beads number
            self.noop = True
        else:
            self.noop = False
            self._b1tob2 = mk_rs_matrix(nbeads1, nbeads2)
            self._b2tob1 = self._b1tob2.T * (float(nbeads1) / float(nbeads2))
            # definition of the scaling also using the open case normal mode matrix transformations
            self._o_b1tob2 = mk_o_rs_matrix(nbeads1, nbeads2)
            self._o_b2tob1 = self._o_b1tob2.T * (float(nbeads1) / float(nbeads2))

    def b1tob2(self, q):
        """Transforms a matrix from one value of beads to another.

        Args:
           q: A matrix with nbeads1 rows, in the bead representation.
        """

        if self.noop:
            # still must return a copy, as the contraction is meant to return new data, not a view
            q_scal = dstrip(q).copy()
        else:
            # this applies to both bead property arrays (e.g. potentials) and bead vector properties (e.g. positions, forces)
            q_scal = np.tensordot(self._b1tob2, q, axes=(1, 0))
            if len(self._open) > 0:
                if len(q_scal.shape) == 2:
                    for io in self._open:
                        q_scal[:, 3 * io] = np.tensordot(
                            self._o_b1tob2, q[:, 3 * io], axes=(1, 0)
                        )
                        q_scal[:, 3 * io + 1] = np.tensordot(
                            self._o_b1tob2, q[:, 3 * io + 1], axes=(1, 0)
                        )
                        q_scal[:, 3 * io + 2] = np.tensordot(
                            self._o_b1tob2, q[:, 3 * io + 2], axes=(1, 0)
                        )
                else:
                    # this applies the open path contraction to EVERYTHING because otherwise we don't know how to handle
                    # the fact that only some beads are open. clearly this is a hack, and in practice the point is that
                    # a "per bead" NM transformation of the potential is not well-defined when different beads have different
                    # NM transformations
                    q_scal = np.tensordot(self._o_b1tob2, q, axes=(1, 0))

        return q_scal

    def b2tob1(self, q):
        """Transforms a matrix from one value of beads to another.

        Args:
           q: A matrix with nbeads2 rows, in the bead representation.
        """

        if self.noop:
            # still must return a copy, as the contraction is meant to return new data, not a view
            q_scal = dstrip(q).copy()
        else:
            # see b1tob2 for an explanation of the rationale for dealing with open path transformations
            q_scal = np.tensordot(self._b2tob1, q, axes=(1, 0))
            if len(self._open) > 0:
                if len(q_scal.shape) == 2:
                    for io in self._open:
                        q_scal[:, 3 * io] = np.tensordot(
                            self._o_b2tob1, q[:, 3 * io], axes=(1, 0)
                        )
                        q_scal[:, 3 * io + 1] = np.tensordot(
                            self._o_b2tob1, q[:, 3 * io + 1], axes=(1, 0)
                        )
                        q_scal[:, 3 * io + 2] = np.tensordot(
                            self._o_b2tob1, q[:, 3 * io + 2], axes=(1, 0)
                        )
                    else:
                        q_scal = np.tensordot(self._o_b2tob1, q, axes=(1, 0))
        return q_scal


class nm_fft(
    object
):  # ! TODO add (matrix-version) of the open path transformation here
    """Uses Fast Fourier transforms to do normal mode transformations.

    Attributes:
       fft: The fast-Fourier transform function to transform between the
          bead and normal mode representations.
       ifft: The inverse fast-Fourier transform function to transform
          between the normal mode and bead representations.
       qdummy: A matrix to hold a copy of the bead positions to transform
          them to the normal mode representation.
       qnmdummy: A matrix to hold a copy of the normal modes to transform
          them to the bead representation.
       nbeads: The number of beads.
       natoms: The number of atoms.
    """

    def __init__(
        self, nbeads, natoms, open_paths=None, n_threads=1, single_precision=False
    ):
        """Initializes nm_trans.

        Args:
           nbeads: The number of beads.
           natoms: The number of atoms.
        """

        self.nbeads = nbeads
        self.natoms = natoms
        self.n_threads = n_threads
        self.single_precision = single_precision
        if open_paths is None:
            open_paths = []
        self._open = open_paths
        # for atoms with open path we still use the matrix transformation
        self._b2o_nm = mk_o_nm_matrix(nbeads)
        self._o_nm2b = self._b2o_nm.T
        try:
            import pyfftw

            info("Import of PyFFTW successful", verbosity.medium)
            self.qdummy = pyfftw.n_byte_align_empty(
                (nbeads, 3 * natoms),
                16,
                "float32" if self.single_precision else "float64",
            )
            self.qnmdummy = pyfftw.n_byte_align_empty(
                (nbeads // 2 + 1, 3 * natoms),
                16,
                "complex64" if self.single_precision else "complex128",
            )

            pyfftw.config.NUM_THREADS = self.n_threads

            self.pyfftw_fw = pyfftw.FFTW(
                self.qdummy,
                self.qnmdummy,
                axes=(0,),
                direction="FFTW_FORWARD",
                threads=self.n_threads,
            )
            self.pyfftw_bw = pyfftw.FFTW(
                self.qnmdummy,
                self.qdummy,
                axes=(0,),
                direction="FFTW_BACKWARD",
                threads=self.n_threads,
            )

            # wrapping these calls into functions has negligible overhead and makes profiling easier
            def call_fft(self):
                self.pyfftw_fw()

            def call_ifft(self):
                self.pyfftw_bw()

            self.fft = lambda: call_fft(self)
            self.ifft = lambda: call_ifft(self)
        except ImportError:  # Uses standard numpy fft library if nothing better
            # is available
            info(
                "Import of PyFFTW unsuccessful, using NumPy library instead. ",
                verbosity.low,
            )
            if self.nbeads > 1:
                warning(
                    "Install PyFFTW to avoid slow normal-modes integration",
                    verbosity.low,
                )

            self.qdummy = np.zeros(
                (nbeads, 3 * natoms),
                dtype="float32" if self.single_precision else "float64",
            )
            self.qnmdummy = np.zeros(
                (nbeads // 2 + 1, 3 * natoms),
                dtype="complex64" if self.single_precision else "complex128",
            )

            def dummy_fft(self):
                self.qnmdummy = np.fft.rfft(self.qdummy, axis=0)

            def dummy_ifft(self):
                self.qdummy = np.fft.irfft(self.qnmdummy, n=self.nbeads, axis=0)

            self.fft = lambda: dummy_fft(self)
            self.ifft = lambda: dummy_ifft(self)

    def b2nm(self, q):
        """Transforms a matrix to the normal mode representation.

        Args:
           q: A matrix with nbeads rows and 3*natoms columns,
              in the bead representation.
        """

        if self.nbeads == 1:
            return q
        self.qdummy[:] = q
        self.fft()
        if self.nbeads == 2:
            return self.qnmdummy.real / np.sqrt(self.nbeads)

        nmodes = self.nbeads // 2

        self.qnmdummy /= np.sqrt(self.nbeads)
        qnm = np.zeros(q.shape)
        qnm[0, :] = self.qnmdummy[0, :].real

        if self.nbeads % 2 == 0:
            self.qnmdummy[1:-1, :] *= np.sqrt(2)
            qnm[1:nmodes, :], qnm[self.nbeads : nmodes : -1, :] = (
                self.qnmdummy[1:-1, :].real,
                self.qnmdummy[1:-1, :].imag,
            )
            qnm[nmodes, :] = self.qnmdummy[nmodes, :].real
        else:
            self.qnmdummy[1:, :] *= np.sqrt(2)
            qnm[1 : nmodes + 1, :], qnm[self.nbeads : nmodes : -1, :] = (
                self.qnmdummy[1:, :].real,
                self.qnmdummy[1:, :].imag,
            )

        for (
            io
        ) in (
            self._open
        ):  # does separately the transformation for the atom that are marked as open paths
            qnm[:, 3 * io] = np.dot(self._b2o_nm, q[:, 3 * io])
            qnm[:, 3 * io + 1] = np.dot(self._b2o_nm, q[:, 3 * io + 1])
            qnm[:, 3 * io + 2] = np.dot(self._b2o_nm, q[:, 3 * io + 2])
        return qnm

    def nm2b(self, qnm):
        """Transforms a matrix to the bead representation.

        Args:
           qnm: A matrix with nbeads rows and 3*natoms columns,
              in the normal mode representation.
        """

        nbeads = self.nbeads
        if nbeads == 1:
            return qnm
        if nbeads == 2:
            self.qnmdummy[:] = qnm
            self.ifft()
            return self.qdummy * np.sqrt(nbeads)

        nmodes = nbeads // 2
        odd = nbeads - 2 * nmodes  # 0 if even, 1 if odd

        isqrt2 = np.sqrt(0.5)
        qnm_complex = self.qnmdummy
        qnm_complex[0, :] = qnm[0, :]
        if not odd:
            qnm_complex[1:-1, :].real = qnm[1:nmodes, :] * isqrt2
            qnm_complex[1:-1, :].imag = qnm[nbeads:nmodes:-1, :] * isqrt2
            qnm_complex[nmodes, :] = qnm[nmodes, :]
        else:
            qnm_complex[1:, :].real = qnm[1 : nmodes + 1, :] * isqrt2
            qnm_complex[1:, :].imag = qnm[nbeads:nmodes:-1, :] * isqrt2

        self.ifft()
        q = self.qdummy * np.sqrt(nbeads)
        for (
            io
        ) in (
            self._open
        ):  # does separately the transformation for the atom that are marked as open paths
            q[:, 3 * io] = np.dot(self._o_nm2b, qnm[:, 3 * io])
            q[:, 3 * io + 1] = np.dot(self._o_nm2b, qnm[:, 3 * io + 1])
            q[:, 3 * io + 2] = np.dot(self._o_nm2b, qnm[:, 3 * io + 2])
        return q
