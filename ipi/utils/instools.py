import numpy as np

from ipi.utils.messages import verbosity, info
from ipi.utils import units
import ipi.utils.mathtools as mt


def banded_hessian(h, im, masses=True, shift=0.001):
    """Given Hessian in the reduced format (h), construct
    the upper band hessian including the RP terms.
    If masses is True returns hessian otherwise it returns dynmat
    shift is value that is added to the diagonal to avoid numerical problems with close to 0 frequencies"""
    nbeads = im.dbeads.nbeads
    natoms = im.dbeads.natoms
    coef = im.coef  # new_disc

    ii = natoms * 3 * nbeads
    ndiag = natoms * 3 + 1  # only upper diagonal form

    href = np.zeros((ndiag, ii))

    # add physical part
    for i in range(nbeads):
        h_aux = h[
            :, i * natoms * 3 : (i + 1) * natoms * 3
        ]  # Peaks one physical hessian
        for j in range(1, ndiag):
            href[j, (ndiag - 1 - j) + i * natoms * 3 : (i + 1) * natoms * 3] = np.diag(
                h_aux, ndiag - 1 - j
            )

    # add spring parts
    if nbeads > 1:
        # Diagonal
        if masses:
            d_corner = im.dbeads.m3[0] * im.omega2
        else:
            d_corner = np.ones(im.dbeads.m3[0].shape) * im.omega2

        d_0 = np.array([[d_corner * 2]]).repeat(im.dbeads.nbeads - 2, axis=0).flatten()
        diag_sp = np.concatenate((d_corner, d_0, d_corner))
        href[-1, :] += diag_sp

        # Non-Diagonal
        d_out = -d_corner
        ndiag_sp = np.array([[d_out]]).repeat(im.dbeads.nbeads - 1, axis=0).flatten()
        href[0, :] = np.concatenate((np.zeros(natoms * 3), ndiag_sp))

    # Add safety shift value
    href[-1, :] += shift

    # ------- new Discretization --------------
    hnew = np.zeros((ndiag, ii))

    # add physical part
    for i in range(nbeads):
        h_aux = (
            h[:, i * natoms * 3 : (i + 1) * natoms * 3] * (coef[i] + coef[i + 1]) / 2
        )  # Peaks one physical hessian
        for j in range(1, ndiag):
            hnew[j, (ndiag - 1 - j) + i * natoms * 3 : (i + 1) * natoms * 3] = np.diag(
                h_aux, ndiag - 1 - j
            )

    if nbeads > 1:
        # Diagonal
        if masses:
            d_corner = im.dbeads.m3[0] * im.omega2
        else:
            d_corner = np.ones(im.dbeads.m3[0].shape) * im.omega2

        d_init = d_corner / coef[1]
        d_fin = d_corner / coef[-2]

        d_mid = d_corner * (1.0 / coef[1] + 1.0 / coef[2])
        for i in range(2, im.dbeads.nbeads - 1):
            d_mid = np.concatenate(
                (d_mid, d_corner * (1.0 / coef[i] + 1.0 / coef[i + 1]))
            )

        diag_sp = np.concatenate((d_init, d_mid, d_fin))
        hnew[-1, :] += diag_sp

        # Non-Diagonal
        d_mid = -d_corner * (1.0 / coef[1])
        for i in range(2, im.dbeads.nbeads):
            d_mid = np.concatenate((d_mid, -d_corner * (1.0 / coef[i])))
        hnew[0, :] = np.concatenate((np.zeros(natoms * 3), d_mid))

    # Add safety shift value
    hnew[-1, :] += shift

    return hnew


def sym_band(A):
    """Returns symmetric banded matrix from just upper banded."""
    u = len(A) - 1
    lu = u
    M = A.shape[1]
    newA = np.empty((u + lu + 1, M))
    newA[: u + 1] = A
    for i in range(1, lu + 1):
        newA[u + i, : M - i] = A[-1 - i, i:]
    return newA


def invmul_banded(A, B, posdef=False):
    """A is in upper banded form
        Solve H.h = -G for Newton - Raphson step, h
    using invmul_banded(H, -G) take step x += h
    to  find minimum or transition state"""

    try:
        from scipy import linalg

        info("Import of scipy successful", verbosity.high)
    except ImportError:
        raise ValueError(" ")

    if posdef:
        return linalg.solveh_banded(A, B)
    else:
        u = len(A) - 1
        lu = u
        newA = sym_band(A)
        # np.set_printoptions(precision=6, suppress=True, threshold=np.nan, linewidth=1000)
        # print linalg.eigvals_banded(A)
        # sys.exit(0)
        return linalg.solve_banded((lu, u), newA, B)


def diag_banded(A, n=2):
    """A is in upper banded form.
    Returns the smallest n eigenvalue and its corresponding eigenvector.
    """
    try:
        from scipy.linalg import eig_banded

        info("Import of scipy successful", verbosity.high)
    except ImportError:
        raise ValueError(" ")

    d = eig_banded(
        A, select="i", select_range=(0, n), eigvals_only=True, check_finite=False
    )

    return d


def red2comp(h, nbeads, natoms, coef=None):
    """Takes the reduced physical hessian (3*natoms*nbeads,3*natoms)
    and construct the 'complete' one (3*natoms*nbeads)^2"""
    info("\n @Instanton: Creating 'complete' physical hessian \n", verbosity.high)

    if coef is None:
        coef = np.ones(nbeads + 1).reshape(-1, 1)
    i = natoms * 3
    ii = nbeads * i
    h0 = np.zeros((ii, ii), float)

    for j in range(nbeads):
        h0[j * i : (j + 1) * i, j * i : (j + 1) * i] = (
            h[:, j * i : (j + 1) * i] * (coef[j] + coef[j + 1]) / 2
        )
    return h0


def get_imvector(h, m3):
    """Compute eigenvector  corresponding to the imaginary mode
    IN     h      = hessian
           m3     = mass vector (dimension = 1 x 3*natoms)
    OUT    imv    = eigenvector corresponding to the imaginary mode
    """
    info("@get_imvector", verbosity.high)
    if h.size != m3.size ** 2:
        raise ValueError(
            "@Get_imvector. Initial hessian size does not match system size."
        )
    m = 1.0 / (m3 ** 0.5)
    mm = np.outer(m, m)
    hm = np.multiply(h, mm)

    # Simmetrize to use linalg.eigh
    hmT = hm.T
    hm = (hmT + hm) / 2.0

    d, w = np.linalg.eigh(hm)
    freq = np.sign(d) * np.absolute(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)

    info(" @GEOP: 1 frequency %4.1f cm^-1" % freq[0], verbosity.low)
    info(" @GEOP: 2 frequency %4.1f cm^-1" % freq[1], verbosity.low)
    info(" @GEOP: 3 frequency %4.1f cm^-1" % freq[2], verbosity.low)
    if freq[0] > -80 and freq[0] < 0:
        raise ValueError(
            " @GEOP: Small negative frequency %4.1f cm^-1" % freq, verbosity.low
        )
    elif freq[0] > 0:
        raise ValueError(
            "@GEOP: The smallest frequency is positive. We aren't in a TS. Please check your hessian"
        )

    info(
        " @get_imvector: We stretch along the mode with freq %f cm^1" % freq[0],
        verbosity.low,
    )

    imv = w[:, 0] * (m3[:] ** 0.5)
    imv = imv / np.linalg.norm(imv)

    return imv.reshape(1, imv.size)


# def print_instanton_geo(prefix, step, nbeads, natoms, names, q, pots, cell, shift):


def print_instanton_geo(
    prefix, step, nbeads, natoms, names, q, f, pots, cell, shift, output_maker
):
    """ Alternative (but very useful) output of the instanton geometry and potential energy """

    outfile = output_maker.get_output(prefix + "_" + str(step) + ".ener", "w")
    print("#Bead    Energy (eV)", file=outfile)
    for i in range(nbeads):
        print(
            (
                str(i)
                + "     "
                + str(units.unit_to_user("energy", "electronvolt", pots[i] - shift))
            ),
            file=outfile,
        )
    outfile.close_stream()

    # print_file("xyz", pos[0], cell, out, title='positions{angstrom}')

    unit = "angstrom"
    unit2 = "atomic_unit"
    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

    outfile = output_maker.get_output(prefix + "_" + str(step) + ".xyz", "w")
    outfile2 = output_maker.get_output(prefix + "_forces_" + str(step) + ".xyz", "w")
    for i in range(nbeads):
        print(natoms, file=outfile)
        print(natoms, file=outfile2)

        print(
            (
                "CELL(abcABC):  %f %f %f %f %f %f cell{atomic_unit}  Traj: positions{%s}   Bead:       %i"
                % (a, b, c, alpha, beta, gamma, unit, i)
            ),
            file=outfile,
        )
        print(
            (
                "CELL(abcABC):  %f %f %f %f %f %f cell{atomic_unit}  Traj: positions{%s}   Bead:       %i"
                % (a, b, c, alpha, beta, gamma, unit2, i)
            ),
            file=outfile2,
        )
        # print >> outfile, ('#Potential (eV):   ' + str(units.unit_to_user('energy', "electronvolt", pots[i] - shift)))

        for j in range(natoms):
            print(
                names[j],
                str(units.unit_to_user("length", unit, q[i, 3 * j])),
                str(units.unit_to_user("length", unit, q[i, 3 * j + 1])),
                str(units.unit_to_user("length", unit, q[i, 3 * j + 2])),
                file=outfile,
            )

        for j in range(natoms):
            print(
                names[j],
                str(units.unit_to_user("force", unit2, f[i, 3 * j])),
                str(units.unit_to_user("force", unit2, f[i, 3 * j + 1])),
                str(units.unit_to_user("force", unit2, f[i, 3 * j + 2])),
                file=outfile2,
            )
    outfile.close_stream()
    outfile2.close_stream()


def print_instanton_hess(prefix, step, hessian, output_maker):
    """ Print physical part of the instanton hessian """
    outfile = output_maker.get_output(prefix + ".hess_" + str(step), "w")
    np.savetxt(outfile, hessian.reshape(1, hessian.size))
    outfile.close_stream()


def ms_pathway(pos, m3):
    """ Compute mass scaled pathway """
    dx = list()
    path = np.zeros(pos.shape[0])
    for i in range(1, pos.shape[0]):
        d_norm = np.linalg.norm((pos[i] - pos[i - 1]) * m3[i] ** 0.5)
        dx.append(d_norm)
        path[i] = np.sum(dx[:i])
    return path
