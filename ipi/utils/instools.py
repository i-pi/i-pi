import numpy as np

from ipi.utils.messages import verbosity, info
from ipi.utils import units
import ipi.utils.mathtools as mt

def banded_hessian(h, im, shift=0.001):
    """Given Hessian in the reduced format (h), construct
    the upper band hessian including the RP terms"""
    nbeads = im.dbeads.nbeads
    natoms = im.dbeads.natoms
    ii = natoms * 3 * nbeads
    ndiag = natoms * 3 + 1  # only upper diagonal form

    # np.set_printoptions(precision=6, suppress=True, threshold=np.nan, linewidth=1000)

    hnew = np.zeros((ndiag, ii))

    # add physical part
    for i in range(nbeads):
        h_aux = h[:, i * natoms * 3:(i + 1) * natoms * 3]  # Peaks one physical hessian
        for j in range(1, ndiag):
            hnew[j, (ndiag - 1 - j) + i * natoms * 3:(i + 1) * natoms * 3] = np.diag(h_aux, ndiag - 1 - j)

    # add spring parts

    if nbeads > 1:
        # Diagonal
        d_corner = im.dbeads.m3[0] * im.omega2
        d_0 = np.array([[d_corner * 2]]).repeat(im.dbeads.nbeads - 2, axis=0).flatten()
        diag_sp = np.concatenate((d_corner, d_0, d_corner))
        hnew[-1, :] += diag_sp

        # Non-Diagonal
        d_out = - d_corner
        ndiag_sp = np.array([[d_out]]).repeat(im.dbeads.nbeads - 1, axis=0).flatten()
        hnew[0, :] = np.concatenate((np.zeros(natoms * 3), ndiag_sp))

    # Add safety shift value
    hnew[-1, :] += shift

    return hnew


def sym_band(A):
    """Return symmetric banded matrix from just upper banded."""
    u = len(A) - 1
    l = u
    M = A.shape[1]
    newA = np.empty((u + l + 1, M))
    newA[:u + 1] = A
    for i in xrange(1, l + 1):
        newA[u + i, :M - i] = A[-1 - i, i:]
    return newA


def invmul_banded(A, B, posdef=False):
    """A is in upper banded form
        Solve H.h = -G for Newton - Raphson step, h
    using invmul_banded(H, -G) take step x += h
    to  find minimum or transition state """

    try:
        from scipy import linalg
        info("Import of scipy successful", verbosity.medium)
    except ImportError:
        raise ValueError(" ")

    if posdef:
        return linalg.solveh_banded(A, B)
    else:
        u = len(A) - 1
        l = u
        newA = sym_band(A)
        # np.set_printoptions(precision=6, suppress=True, threshold=np.nan, linewidth=1000)
        # print linalg.eigvals_banded(A)
        # sys.exit(0)
        return linalg.solve_banded((l, u), newA, B)

def diag_banded(A,n=0):
    """A is in upper banded form.
        Returns the smallest n eigenvalue and its corresponding eigenvector.
        """
    try:
        from scipy.linalg import eig_banded
        info("Import of scipy successful", verbosity.medium)
    except ImportError:
        raise ValueError(" ")

    d , w =  eig_banded(A, select='i', select_range=(0, n))

    return d, w

def red2comp(h, nbeads, natoms):
    """Takes the reduced physical hessian (3*natoms*nbeads,3*natoms)
     and construct the 'complete' one (3*natoms*nbeads)^2 """
    info("\n @Instanton: Creating 'complete' physical hessian \n", verbosity.high)
    i = natoms * 3
    ii = nbeads * i
    h0 = np.zeros((ii, ii), float)

    for j in range(nbeads):
        h0[j * i:(j + 1) * i, j * i:(j + 1) * i] = h[:, j * i:(j + 1) * i]
    return h0

def get_imvector(h, m3):
    """ Compute eigenvector  corresponding to the imaginary mode
            IN     h      = hessian
                   m3     = mass vector (dimension = 1 x 3*natoms)
            OUT    imv    = eigenvector corresponding to the imaginary mode
        """
    info("@get_imvector", verbosity.high)
    if h.size != m3.size ** 2:
        raise ValueError("@Get_imvector. Initial hessian size does not match system size.")
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
        raise ValueError(" @GEOP: Small negative frequency %4.1f cm^-1" % freq, verbosity.low)
    elif freq[0] > 0:
        raise ValueError("@GEOP: The smallest frequency is positive. We aren't in a TS. Please check your hessian")

    info(" @get_imvector: We stretch along the mode with freq %f cm^1" % freq[0], verbosity.low)

    imv = w[:, 0] * (m3[:] ** 0.5)
    imv = imv / np.linalg.norm(imv)

    return imv.reshape(1, imv.size)

def print_instanton_geo(prefix, step, nbeads, natoms, names, q, pots, cell, shift):

    outfile = open(prefix + '_' + str(step) + '.ener', 'w')
    print >> outfile, ('#Bead    Energy (eV)')
    for i in range(nbeads):
        print >> outfile, (str(i) + '     ' + str(units.unit_to_user('energy', "electronvolt", pots[i] - shift)))
    outfile.close()

    # print_file("xyz", pos[0], cell, out, title='positions{angstrom}')

    unit = 'angstrom'
    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

    outfile = open(prefix + '_' + str(step) + '.xyz', 'w')
    for i in range(nbeads):
        print >> outfile, natoms
        # print >> outfile, (('CELL(abcABC): Traj: positions(%s) Bead: %i' %(unit,i) ))
        print >> outfile, ('CELL(abcABC):  %f %f %f %f %f %f cell{atomic_unit}  Traj: positions{%s}   Bead:       %i' % (a, b, c, alpha, beta, gamma, unit, i))
        # print >> outfile, ('#Potential (eV):   ' + str(units.unit_to_user('energy', "electronvolt", pots[i] - shift)))

        for j in range(natoms):
            print >> outfile, names[j], \
                str(units.unit_to_user('length', unit, q[i, 3 * j])), \
                str(units.unit_to_user('length', unit, q[i, 3 * j + 1])), \
                str(units.unit_to_user('length', unit, q[i, 3 * j + 2]))

    outfile.close()


def print_instanton_hess(prefix, step, hessian):

    np.set_printoptions(precision=7, suppress=True, threshold=np.nan, linewidth=3000)
    outfile = open(prefix + '.hess_' + str(step), 'w')
    np.savetxt(outfile, hessian.reshape(1, hessian.size))
    outfile.close()
