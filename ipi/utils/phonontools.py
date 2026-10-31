import numpy as np

from ipi.utils.array_backend import xp


def apply_asr(asr, dm, beads, return_trans_matrix=False):
    """
    Removes the translations and/or rotations depending on the asr mode.
    """

    if asr == "none":
        return dm

    ism = 1 / xp.sqrt(beads.m3[-1])
    m = beads.m

    # Computes the centre of mass.
    com = (beads.q.reshape((beads.natoms, 3)).T @ m) / m.sum()
    qminuscom = beads.q.reshape((beads.natoms, 3)) - com
    # Computes the moment of inertia tensor.
    moi = xp.zeros((3, 3))
    for k in range(beads.natoms):
        moi -= (
            xp.linalg.cross(qminuscom[k], xp.eye(3))
            @ xp.linalg.cross(qminuscom[k], xp.eye(3))
            * m[k]
        )

    # eigh gives orthonormal eigenvectors, also for equal moments of inertia
    U = xp.linalg.eigh(moi)[1]
    R = qminuscom @ U

    if asr == "crystal":
        D = xp.zeros((3, 3 * beads.natoms))

        # Computes the vectors along rotations.
        D[0] = xp.asarray(np.tile([1, 0, 0], beads.natoms)) / ism
        D[1] = xp.asarray(np.tile([0, 1, 0], beads.natoms)) / ism
        D[2] = xp.asarray(np.tile([0, 0, 1], beads.natoms)) / ism

        # Computes unit vecs.
        for k in range(3):
            D[k] = D[k] / xp.linalg.vector_norm(D[k])

    elif asr == "poly":
        D = xp.zeros((6, 3 * beads.natoms))

        # Computes the vectors along rotations.
        D[0] = xp.asarray(np.tile([1, 0, 0], beads.natoms)) / ism
        D[1] = xp.asarray(np.tile([0, 1, 0], beads.natoms)) / ism
        D[2] = xp.asarray(np.tile([0, 0, 1], beads.natoms)) / ism

        # element (iatom, idof) goes to D[:, 3 * iatom + idof]
        D[3] = (R[:, 1:2] * U[:, 2] - R[:, 2:3] * U[:, 1]).reshape(-1) / ism
        D[4] = (R[:, 2:3] * U[:, 0] - R[:, 0:1] * U[:, 2]).reshape(-1) / ism
        D[5] = (R[:, 0:1] * U[:, 1] - R[:, 1:2] * U[:, 0]).reshape(-1) / ism

        # Computes unit vecs.
        for k in range(6):
            D[k] = D[k] / xp.linalg.vector_norm(D[k])

    # Computes the transformation matrix.
    transfmatrix = xp.eye(3 * beads.natoms) - D.T @ D

    if return_trans_matrix:
        return transfmatrix
    else:
        return transfmatrix.T @ (dm @ transfmatrix)
