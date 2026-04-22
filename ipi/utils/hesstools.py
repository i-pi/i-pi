"""Functions that deal with hessian transformations and calculations"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import os

import numpy as np

from ipi.utils.array_backend import xp, xp_size, to_numpy
from ipi.utils.depend import dstrip
from ipi.utils.messages import verbosity, info


def get_dynmat(h, m3, nbeads=1):
    """Computes the dynamical matrix.
    If nbeads > 1 the reduced form with shape = (3*natoms, 3*natoms*nbeads) is expected
    """

    # Check dimensions

    if h.shape != (m3.shape[1], m3.shape[1] * nbeads):
        print(h.shape, m3.shape)
        raise ValueError(
            "@get_dynmat: The provided hessian hasn't the proper dimension (3*natoms, 3*natoms*nbeads) "
        )

    m3 = dstrip(m3)
    h = dstrip(h)
    ism = xp.reshape(m3, (1, -1)) ** (-0.5)
    ismT = xp.reshape(m3[0], (-1, 1)) ** (-0.5)

    dynmat = ismT * (h * ism)
    return dynmat


def clean_hessian(h, q, natoms, nbeads, m, m3, asr, mofi=False):
    """
    Removes the translations and rotations modes.
    IN  h      = hessian (3*natoms*nbeads, 3*natoms*nbeads)
        q      = positions
        natoms = number of atoms
        nbeads = number of beads
        m      = mass vector, one value for each atom
        m3     = mass vector, one value for each degree of freedom
        mofi   = An optional boolean which decides whether the det(M_of_I)
                 is returned or not
    OUT d      = non zero eigenvalues of the dynmatrix
        w      = eigenvectors without external modes

    #Adapted from ipi/engine/motion/phonons.py apply_asr"""

    info(" @clean_hessian: asr = %s " % asr, verbosity.medium)
    # Inputs may arrive wrapped in depend_arrays; peel them once.
    m = dstrip(m)
    m3 = dstrip(m3)
    q = dstrip(q)
    h = dstrip(h)
    # Set some useful things
    ii = natoms * nbeads
    mm = xp.tile(m, (nbeads,))
    ism = xp.reshape(m3, (ii * 3, 1)) ** (-0.5)
    dynmat = ism.T * (h * ism)

    if asr == "none" or asr is None:
        hm = dynmat
    else:
        # Computes the centre of mass.
        qr = xp.reshape(q, (ii, 3))
        com = (xp.matrix_transpose(qr) @ mm) / xp.sum(mm)
        qminuscom = qr - com
        ism = xp.reshape(ism, (-1,))

        if asr == "poly":
            # Computes the moment of inertia tensor.
            moi = xp.zeros((3, 3))
            eye3 = xp.eye(3)
            for k in range(ii):
                cx = xp.linalg.cross(qminuscom[k], eye3)
                moi = moi - (cx @ cx) * mm[k]

            I, U = xp.linalg.eigh(moi)
            R = qminuscom @ U
            D = xp.zeros((6, 3 * ii))

            # Computes the vectors along translations and rotations.
            # Translations
            ex = xp.asarray(np.tile([1.0, 0.0, 0.0], ii))
            ey = xp.asarray(np.tile([0.0, 1.0, 0.0], ii))
            ez = xp.asarray(np.tile([0.0, 0.0, 1.0], ii))
            D[0] = ex / ism
            D[1] = ey / ism
            D[2] = ez / ism
            # Rotations
            for i in range(3 * ii):
                iatom = i // 3
                idof = i % 3
                D[3, i] = (R[iatom, 1] * U[idof, 2] - R[iatom, 2] * U[idof, 1]) / ism[i]
                D[4, i] = (R[iatom, 2] * U[idof, 0] - R[iatom, 0] * U[idof, 2]) / ism[i]
                D[5, i] = (R[iatom, 0] * U[idof, 1] - R[iatom, 1] * U[idof, 0]) / ism[i]

            for k in range(6):
                D[k] = D[k] / xp.linalg.vector_norm(D[k])
            # Computes the transformation matrix.
            transfmatrix = xp.eye(3 * ii) - (D.T @ D)
            hm = transfmatrix.T @ (dynmat @ transfmatrix)

        elif asr == "crystal":
            # Computes the vectors along translations.
            D = xp.zeros((3, 3 * ii))
            ex = xp.asarray(np.tile([1.0, 0.0, 0.0], ii))
            ey = xp.asarray(np.tile([0.0, 1.0, 0.0], ii))
            ez = xp.asarray(np.tile([0.0, 0.0, 1.0], ii))
            D[0] = ex / ism
            D[1] = ey / ism
            D[2] = ez / ism

            for k in range(3):
                D[k] = D[k] / xp.linalg.vector_norm(D[k])
            # Computes the transformation matrix.
            transfmatrix = xp.eye(3 * ii) - (D.T @ D)
            hm = transfmatrix.T @ (dynmat @ transfmatrix)

    # Symmetrize to use linalg.eigh
    hm = (hm.T + hm) / 2.0

    d, w = xp.linalg.eigh(hm)

    # convert to cm^-1
    dd = xp.sign(d) * xp.abs(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)

    # Zeros
    cut0 = 0.01  # Note that dd[] units are cm^1
    nzero = dd[xp.abs(dd) < cut0]
    n_zero = int(xp_size(nzero))

    if asr == "poly" and n_zero != 6:
        info(" @GEOP: Warning, we have %d 'zero' frequencies" % n_zero, verbosity.low)

    if asr == "crystal" and n_zero != 3:
        info(" @GEOP: Warning, we have %d 'zero' frequencies" % n_zero, verbosity.low)

    # Negatives
    cutNeg = -4  # Note that dd[] units are cm^1
    nneg = dd[dd < cutNeg]
    n_neg = int(xp_size(nneg))
    info(" @Clean hessian: We have %d 'neg' frequencies " % n_neg, verbosity.medium)

    # Now eliminate external degrees of freedom from the dynmatrix

    if n_zero > 0:
        if float(xp.linalg.vector_norm(nzero)) > cut0:
            info(
                " Warning @Clean hessian: We have deleted %d 'zero' frequencies "
                % n_zero,
                verbosity.high,
            )
            info(
                " but the norm is greater than 0.01 cm^-1.  This should not happen.",
                verbosity.high,
            )

        keep = list(range(n_neg)) + list(range(n_neg + n_zero, int(xp_size(d))))
        d = d[keep]
        w = w[:, keep]

    if mofi:
        if asr == "poly":
            return d, w, xp.prod(I)
        else:
            return d, w, 1.0
    else:
        return d, w


def get_hessian(
    gm, x0, natoms, nbeads=1, fixatoms_dof=[], d=0.001, new_disc=False, friction=False
):
    """Compute the physical hessian given a function to evaluate energy and forces (gm).
    The intermediate steps are written as a temporary files so the full hessian calculations is only ONE step.

    IN     gm       = gradient mapper
           x0       = position vector
           natoms   = number of atoms
           nbeads   = number of beads
           fixatoms_dof = indexes of fixed degrees of freedom
           d        = displacement

    OUT    h        = physical hessian  (3*natoms-len(fixatoms_dof)  , nbeads*( 3*natoms-len(fixatoms_dof) ))
    """

    info(" @get_hessian: Computing hessian", verbosity.low)
    ii = natoms * 3
    activedof = np.delete(np.arange(ii), fixatoms_dof)
    ncalc = ii - len(fixatoms_dof)
    if xp_size(x0) != natoms * 3 * nbeads:
        raise ValueError(
            "The position vector is not consistent with the number of atoms/beads."
        )

    h = xp.zeros((ii, ii * nbeads))
    if friction:
        eta_h = xp.zeros((nbeads, ii, ii, ii))

    # Check if there is a temporary file:
    i0 = -1

    for i in range(ii, -1, -1):
        try:
            b = np.loadtxt("hessian_" + str(i) + ".tmp")
        except IOError:
            pass
        else:
            h[:, :] = b[:, :]
            i0 = i
            print(("We have found a temporary file ( hessian_" + str(i) + ".tmp). "))
            if (
                b.shape == h.shape
            ):  # Check that the last temporary file was properly written
                break
            else:
                continue
    if friction:
        for i in range(ii, -1, -1):
            try:
                b = np.loadtxt("hessianEta_" + str(i) + ".tmp")
            except IOError:
                pass
            else:
                eta_h[:] = b.reshape((nbeads, ii, ii, ii))
                i0 = i
                print(
                    (
                        "We have found a temporary file ( hessianEta_"
                        + str(i)
                        + ".tmp). "
                    )
                )
                if (
                    b.shape == eta_h.shape
                ):  # Check that the last temporary file was properly written
                    break
                else:
                    continue

    # Start calculation:
    for j in range(i0 + 1, ii):
        if j in fixatoms_dof:
            continue
        else:
            ndone = len(activedof[activedof < j])
            info(
                " @get_hessian: Computing hessian: %d of %d" % (ndone + 1, ncalc),
                verbosity.low,
            )
            x = xp.asarray(x0, copy=True)

            # PLUS
            x[:, j] = x0[:, j] + d
            e, f1 = gm(x, new_disc)
            if friction:
                eta1 = gm.eta

            # Minus
            x[:, j] = x0[:, j] - d
            e, f2 = gm(x, new_disc)
            if friction:
                eta2 = gm.eta

            # COMBINE
            g = (f1 - f2) / (2 * d)
            h[j, :] = g.flatten()
            f = open("hessian_" + str(j) + ".tmp", "w")
            np.savetxt(f, to_numpy(h))
            f.close()

            if friction:
                eta_h[:, :, :, j] = (eta1 - eta2) / (2 * d)
                f = open("hessianEta_" + str(j) + ".tmp", "w")
                np.savetxt(f, to_numpy(eta_h).flatten())
                f.close()

    u, g = gm(x0)  # Keep the mapper updated

    for i in range(ii):
        try:
            os.remove("hessianEta_" + str(i) + ".tmp")
        except OSError:
            pass

        try:
            os.remove("hessian_" + str(i) + ".tmp")
        except OSError:
            pass

    if not friction:
        return h
    else:
        return h, eta_h
