import numpy as np

from ipi.utils.array_backend import xp


def apply_asr(asr, dm, beads, return_trans_matrix=False):
    """
    Removes the translations and/or rotations depending on the asr mode.
    """

    if asr == "none":
        return dm

    nat = beads.natoms
    ism = 1 / xp.sqrt(beads.m3[-1])
    m = beads.m

    # Computes the centre of mass.
    qr = xp.reshape(beads.q, (nat, 3))
    com = (xp.matrix_transpose(qr) @ m) / xp.sum(m)
    qminuscom = qr - com
    # Computes the moment of inertia tensor.
    moi = xp.zeros((3, 3))
    eye3 = xp.eye(3)
    for k in range(nat):
        cx = xp.linalg.cross(qminuscom[k], eye3)
        moi = moi - (cx @ cx) * m[k]

    U = xp.linalg.eig(moi)[1]
    R = qminuscom @ U

    ex = xp.asarray(np.tile([1.0, 0.0, 0.0], nat))
    ey = xp.asarray(np.tile([0.0, 1.0, 0.0], nat))
    ez = xp.asarray(np.tile([0.0, 0.0, 1.0], nat))

    if asr == "crystal":
        D = xp.zeros((3, 3 * nat))

        D[0] = ex / ism
        D[1] = ey / ism
        D[2] = ez / ism

        for k in range(3):
            D[k] = D[k] / xp.linalg.vector_norm(D[k])

    elif asr == "poly":
        D = xp.zeros((6, 3 * nat))

        D[0] = ex / ism
        D[1] = ey / ism
        D[2] = ez / ism

        for i in range(3 * nat):
            iatom = i // 3
            idof = i % 3
            D[3, i] = (R[iatom, 1] * U[idof, 2] - R[iatom, 2] * U[idof, 1]) / ism[i]
            D[4, i] = (R[iatom, 2] * U[idof, 0] - R[iatom, 0] * U[idof, 2]) / ism[i]
            D[5, i] = (R[iatom, 0] * U[idof, 1] - R[iatom, 1] * U[idof, 0]) / ism[i]

        for k in range(6):
            D[k] = D[k] / xp.linalg.vector_norm(D[k])

    # Computes the transformation matrix.
    transfmatrix = xp.eye(3 * nat) - (D.T @ D)

    if return_trans_matrix:
        return transfmatrix
    else:
        return transfmatrix.T @ (dm @ transfmatrix)
