from ipi.utils.array_backend import xp


def external_modes(q, m, asr):
    """Returns the orthonormal mass-scaled displacements along the
    translations of points of masses `m` at positions `q` (one row per
    point), and also along their rotations if `asr` is "poly", as the rows of
    an array, and the principal moments of inertia ("poly") or None."""

    if asr not in ("crystal", "poly"):
        raise ValueError(f"Unknown acoustic sum rule '{asr}'")

    sm = xp.sqrt(m)[:, None]
    eye = xp.eye(3)
    # translations
    modes = [(sm * eye[k]).reshape(-1) for k in range(3)]
    inertia = None
    if asr == "poly":
        r = q - (q.T @ m) / m.sum()
        # the moment of inertia tensor, sum_k m_k (|r_k|^2 I - r_k r_k^T)
        moi = eye * (m @ (r * r).sum(axis=1)) - (r.T * m) @ r
        # eigh gives orthonormal eigenvectors, also for equal moments of inertia
        inertia, U = xp.linalg.eigh(moi)
        R = r @ U
        # rotations around the principal axes
        modes += [
            ((R[:, 1:2] * U[:, 2] - R[:, 2:3] * U[:, 1]) * sm).reshape(-1),
            ((R[:, 2:3] * U[:, 0] - R[:, 0:1] * U[:, 2]) * sm).reshape(-1),
            ((R[:, 0:1] * U[:, 1] - R[:, 1:2] * U[:, 0]) * sm).reshape(-1),
        ]
    D = xp.stack(modes)
    return D / xp.linalg.vector_norm(D, axis=1)[:, None], inertia


def apply_asr(asr, dm, beads, return_trans_matrix=False):
    """
    Removes the translations and/or rotations depending on the asr mode.
    """

    if asr == "none":
        return dm

    D = external_modes(beads.q.reshape((beads.natoms, 3)), beads.m, asr)[0]
    # projector on the space orthogonal to the external modes
    transfmatrix = xp.eye(3 * beads.natoms) - D.T @ D

    if return_trans_matrix:
        return transfmatrix
    else:
        return transfmatrix.T @ (dm @ transfmatrix)
