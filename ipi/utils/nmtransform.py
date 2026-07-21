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
