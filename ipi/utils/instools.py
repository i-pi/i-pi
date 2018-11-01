import numpy as np

from ipi.utils.messages import verbosity, info
from ipi.utils import units
import ipi.utils.mathtools as mt
import os.path


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


def red2comp(h, nbeads, natoms):
    """Takes the reduced physical hessian and construct the 'complete' one (all 0 included) """
    info("\n @Instanton: Creating 'complete' physical hessian \n", verbosity.high)
    i = natoms * 3
    ii = nbeads * i
    h0 = np.zeros((ii, ii), float)

    for j in range(nbeads):
        h0[j * i:(j + 1) * i, j * i:(j + 1) * i] = h[:, j * i:(j + 1) * i]
    return h0


def get_hessian(h, gm, x0, d=0.0005):
    """Compute the physical hessian
       IN     h       = physical hessian
              gm      = gradient mapper
              x0      = position vector

       OUT    h       = physical hessian
        """
    # TODO What about the case you have numerical gradients?

    info(" @Instanton: Computing hessian", verbosity.low)
    ii = gm.dbeads.natoms * 3
    h[:] = np.zeros((h.shape), float)

    # Ask Michele about transfer force here I
    # ddbeads = gm.dbeads.copy()
    # ddcell = gm.dcell.copy()
    # ddforces = gm.dforces.copy(ddbeads, ddcell)

    i0 = -1
    # Check if there is a temporal file:
    for i in range(ii, -1, -1):
        try:
            b = np.loadtxt('hessian_' + str(i) + '.tmp')
        except IOError:
            pass
        else:
            h[:, :] = b[:, :]
            i0 = i
            print('We have found a temporary file ( hessian_' + str(i) + '.tmp). ')
            if b.shape == h.shape:  # Check that the last temporal file was properly written
                break
            else:
                continue

    for j in range(i0 + 1, ii):
        info(" @Instanton: Computing hessian: %d of %d" % ((j + 1), ii), verbosity.low)
        x = x0.copy()

        x[:, j] = x0[:, j] + d
        e, f1 = gm(x)
        x[:, j] = x0[:, j] - d
        e, f2 = gm(x)
        g = (f1 - f2) / (2 * d)

        h[j, :] = g.flatten()

        f = open('hessian_' + str(j) + '.tmp', 'w')
        #print >> f, 'STEP %i' % j
        np.savetxt(f, h)
        f.close()

    u, g = gm(x0)  # Keep the mapper updated

    # Ask Michele about transfer force here II
    # gm.dbeads.q = ddbeads.q
    # gm.dforces.transfer_forces(ddforces)

    for i in range(ii):
        try:
            os.remove('hessian_' + str(i) + '.tmp')
        except OSError:
            pass


def clean_hessian(h, q, natoms, nbeads, m, m3, asr, mofi=False):
    """
        Removes the translations and rotations modes.
        IN  h      = hessian
            q      = positions
            natoms = number of atoms
            nbeads = number of beads
            m      = mass vector, one value for each atom
            m3     = mass vector, one value for each degree of freedom
            mofi   = An optional boolean which decides whether the det(M_of_I)
                     is returned or not
        OUT d      = non zero eigenvalues of the dynmatrix
            w      = eigenvectors without external modes

        #Adapted from ipi/engine/motion/phonons.py apply_asr    """

    info(" @clean_hessian", verbosity.high)
    # Set some useful things
    ii = natoms * nbeads
    mm = np.zeros((nbeads, natoms))
    for i in range(nbeads):
        mm[i] = m
    mm = mm.reshape(ii)
    ism = m3.reshape(ii * 3) ** (-0.5)
    ismm = np.outer(ism, ism)
    dynmat = np.multiply(h, ismm)

    if asr == 'none':
        hm = dynmat
    else:
        # Computes the centre of mass.
        com = np.dot(np.transpose(q.reshape((ii, 3))), mm) / mm.sum()
        qminuscom = q.reshape((ii, 3)) - com

        if asr == 'poly':
            # Computes the moment of inertia tensor.
            moi = np.zeros((3, 3), float)
            for k in range(ii):
                moi -= np.dot(np.cross(qminuscom[k], np.identity(3)), np.cross(qminuscom[k], np.identity(3))) * mm[k]

            I, U = (np.linalg.eig(moi))
            R = np.dot(qminuscom, U)
            D = np.zeros((6, 3 * ii), float)

            # Computes the vectors along translations and rotations.
            # Translations
            D[0] = np.tile([1, 0, 0], ii) / ism
            D[1] = np.tile([0, 1, 0], ii) / ism
            D[2] = np.tile([0, 0, 1], ii) / ism
            # Rotations
            for i in range(3 * ii):
                iatom = i / 3
                idof = np.mod(i, 3)
                D[3, i] = (R[iatom, 1] * U[idof, 2] - R[iatom, 2] * U[idof, 1]) / ism[i]
                D[4, i] = (R[iatom, 2] * U[idof, 0] - R[iatom, 0] * U[idof, 2]) / ism[i]
                D[5, i] = (R[iatom, 0] * U[idof, 1] - R[iatom, 1] * U[idof, 0]) / ism[i]

            for k in range(6):
                D[k] = D[k] / np.linalg.norm(D[k])
            # Computes the transformation matrix.
            transfmatrix = np.eye(3 * ii) - np.dot(D.T, D)
            hm = np.dot(transfmatrix.T, np.dot(dynmat, transfmatrix))

        elif asr == 'crystal':
            # Computes the vectors along translations.
            # Translations
            D = np.zeros((3, 3 * ii), float)
            D[0] = np.tile([1, 0, 0], ii) / ism
            D[1] = np.tile([0, 1, 0], ii) / ism
            D[2] = np.tile([0, 0, 1], ii) / ism

            for k in range(3):
                D[k] = D[k] / np.linalg.norm(D[k])
            # Computes the transformation matrix.
            transfmatrix = np.eye(3 * ii) - np.dot(D.T, D)
            hm = np.dot(transfmatrix.T, np.dot(dynmat, transfmatrix))

    # Simmetrize to use linalg.eigh
    hmT = hm.T
    hm = (hmT + hm) / 2.0

    d, w = np.linalg.eigh(hm)

    # Count
    dd = np.sign(d) * np.absolute(d) ** 0.5 / (2 * np.pi * 3e10 * 2.4188843e-17)  # convert to cm^-1

    # Zeros
    cut0 = 0.01  # Note that dd[] units are cm^1
    condition = np.abs(dd) < cut0
    nzero = np.extract(condition, dd)

    if asr == 'poly' and nzero.size != 6:
        info(" @GEOP: Warning, we have %d 'zero' frequencies" % nzero.size, verbosity.low)

    if asr == 'crystal' and nzero.size != 3:
        info(" @GEOP: Warning, we have %d 'zero' frequencies" % nzero.size, verbosity.low)

    # Negatives
    cutNeg = -4  # Note that dd[] units are cm^1
    condition = dd < cutNeg
    nneg = np.extract(condition, dd)
    info(" @Clean hessian: We have %d 'neg' frequencies " % (nneg.size), verbosity.medium)

    # Now eliminate external degrees of freedom from the dynmatrix

    if nzero.size > 0:
        if np.linalg.norm(nzero) > cut0:
            info(" Warning @Clean hessian: We have deleted %d 'zero' frequencies " % (nzero.size), verbosity.high)
            info(" but the norm is greater than 0.01 cm^-1.  This should not happen.", verbosity.high)

        d = np.delete(d, range(nneg.size, nneg.size + nzero.size))
        w = np.delete(w, range(nneg.size, nneg.size + nzero.size), axis=1)

    if mofi:
        if asr == 'poly':
            return d, w, np.prod(I)
        else:
            return d, w, 1.0
    else:
        return d, w


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

    return imv


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
