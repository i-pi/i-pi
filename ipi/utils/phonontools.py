import numpy as np


def apply_asr(asr, dm, beads, return_trans_matrix=False):
    """
    Removes the translations and/or rotations depending on the asr mode.
    """

    if asr == "none":
        return dm

    ism = 1 / np.sqrt(beads.m3[-1])
    m = beads.m

    # Computes the centre of mass.
    com = np.dot(np.transpose(beads.q.reshape((beads.natoms, 3))), m) / m.sum()
    qminuscom = beads.q.reshape((beads.natoms, 3)) - com
    # Computes the moment of inertia tensor.
    moi = np.zeros((3, 3), float)
    for k in range(beads.natoms):
        moi -= (
            np.dot(
                np.cross(qminuscom[k], np.identity(3)),
                np.cross(qminuscom[k], np.identity(3)),
            )
            * m[k]
        )

    U = (np.linalg.eig(moi))[1]
    R = np.dot(qminuscom, U)

    if asr == "crystal":
        D = np.zeros((3, 3 * beads.natoms), float)

        # Computes the vectors along rotations.
        D[0] = np.tile([1, 0, 0], beads.natoms) / ism
        D[1] = np.tile([0, 1, 0], beads.natoms) / ism
        D[2] = np.tile([0, 0, 1], beads.natoms) / ism

        # Computes unit vecs.
        for k in range(3):
            D[k] = D[k] / np.linalg.norm(D[k])

    elif asr == "poly":
        D = np.zeros((6, 3 * beads.natoms), float)

        # Computes the vectors along rotations.
        D[0] = np.tile([1, 0, 0], beads.natoms) / ism
        D[1] = np.tile([0, 1, 0], beads.natoms) / ism
        D[2] = np.tile([0, 0, 1], beads.natoms) / ism

        for i in range(3 * beads.natoms):
            iatom = i // 3
            idof = np.mod(i, 3)
            D[3, i] = (R[iatom, 1] * U[idof, 2] - R[iatom, 2] * U[idof, 1]) / ism[i]
            D[4, i] = (R[iatom, 2] * U[idof, 0] - R[iatom, 0] * U[idof, 2]) / ism[i]
            D[5, i] = (R[iatom, 0] * U[idof, 1] - R[iatom, 1] * U[idof, 0]) / ism[i]

        # Computes unit vecs.
        for k in range(6):
            D[k] = D[k] / np.linalg.norm(D[k])

    # Computes the transformation matrix.
    transfmatrix = np.eye(3 * beads.natoms) - np.dot(D.T, D)

    if return_trans_matrix:
        return transfmatrix
    else:
        return np.dot(transfmatrix.T, np.dot(dm, transfmatrix))
