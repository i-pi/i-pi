"""Functions that deal with hessian transformations and calculations"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
from ipi.utils.messages import verbosity, info
from ipi.utils.array_backend import xp, xp_size, to_numpy
import os


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

    ism = m3.reshape((1, -1)) ** (-0.5)
    ismT = m3[0].reshape((-1, 1)) ** (-0.5)

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
    # Set some useful things
    ii = natoms * nbeads
    mm = xp.zeros((nbeads, natoms))
    for i in range(nbeads):
        mm[i] = m
    mm = mm.reshape(ii)
    ism = m3.reshape((ii * 3, 1)) ** (-0.5)
    dynmat = ism.T * (h * ism)
    # ismm = np.outer(ism, ism)
    # dynmat = np.multiply(h, ismm)

    if asr == "none" or asr is None:
        hm = dynmat
    else:
        # Computes the centre of mass.
        com = (q.reshape((ii, 3)).T @ mm) / mm.sum()
        qminuscom = q.reshape((ii, 3)) - com
        ism = ism.reshape(-1)

        if asr == "poly":
            # Computes the moment of inertia tensor.
            moi = xp.zeros((3, 3))
            for k in range(ii):
                moi -= (
                    xp.linalg.cross(qminuscom[k], xp.eye(3))
                    @ xp.linalg.cross(qminuscom[k], xp.eye(3))
                    * mm[k]
                )

            # eigh gives orthonormal eigenvectors, also for equal moments of inertia
            I, U = xp.linalg.eigh(moi)
            R = qminuscom @ U
            D = xp.zeros((6, 3 * ii))

            # Computes the vectors along translations and rotations.
            # Translations
            D[0] = xp.asarray(np.tile([1, 0, 0], ii)) / ism
            D[1] = xp.asarray(np.tile([0, 1, 0], ii)) / ism
            D[2] = xp.asarray(np.tile([0, 0, 1], ii)) / ism
            # Rotations, element (iatom, idof) goes to D[:, 3 * iatom + idof]
            D[3] = (R[:, 1:2] * U[:, 2] - R[:, 2:3] * U[:, 1]).reshape(-1) / ism
            D[4] = (R[:, 2:3] * U[:, 0] - R[:, 0:1] * U[:, 2]).reshape(-1) / ism
            D[5] = (R[:, 0:1] * U[:, 1] - R[:, 1:2] * U[:, 0]).reshape(-1) / ism

            for k in range(6):
                D[k] = D[k] / xp.linalg.vector_norm(D[k])
            # Computes the transformation matrix.
            transfmatrix = xp.eye(3 * ii) - D.T @ D
            hm = transfmatrix.T @ (dynmat @ transfmatrix)

        elif asr == "crystal":
            # Computes the vectors along translations.
            # Translations
            D = xp.zeros((3, 3 * ii))
            D[0] = xp.asarray(np.tile([1, 0, 0], ii)) / ism
            D[1] = xp.asarray(np.tile([0, 1, 0], ii)) / ism
            D[2] = xp.asarray(np.tile([0, 0, 1], ii)) / ism

            for k in range(3):
                D[k] = D[k] / xp.linalg.vector_norm(D[k])
            # Computes the transformation matrix.
            transfmatrix = xp.eye(3 * ii) - D.T @ D
            hm = transfmatrix.T @ (dynmat @ transfmatrix)

    # Symmetrize to use linalg.eigh
    hmT = hm.T
    hm = (hmT + hm) / 2.0

    d, w = xp.linalg.eigh(hm)

    # Count, on the host
    d_host = to_numpy(d)
    dd = (
        np.sign(d_host)
        * np.absolute(d_host) ** 0.5
        / (2 * np.pi * 3e10 * 2.4188843e-17)
    )  # convert to cm^-1

    # Zeros
    cut0 = 0.01  # Note that dd[] units are cm^1
    condition = np.abs(dd) < cut0
    nzero = np.extract(condition, dd)

    if asr == "poly" and nzero.size != 6:
        info(
            " @GEOP: Warning, we have %d 'zero' frequencies" % nzero.size, verbosity.low
        )

    if asr == "crystal" and nzero.size != 3:
        info(
            " @GEOP: Warning, we have %d 'zero' frequencies" % nzero.size, verbosity.low
        )

    # Negatives
    cutNeg = -4  # Note that dd[] units are cm^1
    condition = dd < cutNeg
    nneg = np.extract(condition, dd)
    info(
        " @Clean hessian: We have %d 'neg' frequencies " % (nneg.size), verbosity.medium
    )

    # Now eliminate external degrees of freedom from the dynmatrix

    if nzero.size > 0:
        if np.linalg.norm(nzero) > cut0:
            info(
                " Warning @Clean hessian: We have deleted %d 'zero' frequencies "
                % (nzero.size),
                verbosity.high,
            )
            info(
                " but the norm is greater than 0.01 cm^-1.  This should not happen.",
                verbosity.high,
            )

        # removes the eigenpairs from nneg.size to nneg.size + nzero.size
        d = xp.concat((d[: nneg.size], d[nneg.size + nzero.size :]))
        w = xp.concat((w[:, : nneg.size], w[:, nneg.size + nzero.size :]), axis=1)

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

    # the hessian is assembled on the host, as it is saved to file after each row
    h = np.zeros((ii, ii * nbeads), float)
    if friction:
        eta_h = np.zeros((nbeads, ii, ii, ii), float)

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
            h[j, :] = to_numpy(g).reshape(-1)
            f = open("hessian_" + str(j) + ".tmp", "w")
            np.savetxt(f, h)
            f.close()

            if friction:
                eta_h[:, :, :, j] = (eta1 - eta2) / (2 * d)
                f = open("hessianEta_" + str(j) + ".tmp", "w")
                np.savetxt(f, eta_h.flatten())
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

    # the friction tensors (gm.eta) are kept on the host
    if not friction:
        return xp.asarray(h)
    else:
        return xp.asarray(h), eta_h
