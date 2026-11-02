"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.
"""

__all__ = ["DynMatrixMover"]

import numpy as np


from ipi.engine.motion import Motion
from ipi.engine.motion.motion import zeros_or_check
from ipi.utils.messages import verbosity, info
from ipi.utils.array_backend import xp, xp_size, to_numpy
from ipi.utils.phonontools import apply_asr


class DynMatrixMover(Motion):
    """Dynamic matrix calculation routine by finite difference."""

    jax_support = False

    def __init__(
        self,
        fixcom=False,
        fixatoms_dof=None,
        mode="fd",
        energy_shift=0.0,
        pos_shift=0.001,
        output_shift=0.000,
        dynmat=np.zeros(0, float),
        refdynmat=np.zeros(0, float),
        prefix="",
        asr="none",
    ):
        """Initialises DynMatrixMover.
        Args:
        fixcom  : An optional boolean which decides whether the centre of mass
                  motion will be constrained or not. Defaults to False.
        dynmatrix : A 3Nx3N array that stores the dynamic matrix.
        refdynmatrix : A 3Nx3N array that stores the refined dynamic matrix.
        """

        super(DynMatrixMover, self).__init__(fixcom=fixcom, fixatoms_dof=fixatoms_dof)

        # Finite difference option.
        self.mode = mode
        if self.mode == "fd":
            self.phcalc = FDPhononCalculator()
        elif self.mode == "nmfd":
            self.phcalc = NMFDPhononCalculator()
        elif self.mode == "enmfd":
            self.phcalc = ENMFDPhononCalculator()

        self.deltaw = output_shift
        self.deltax = pos_shift
        self.deltae = energy_shift
        self.dynmatrix = dynmat
        self.refdynmatrix = refdynmat
        self.U = None
        self.V = None
        self.prefix = prefix
        self.asr = asr

        if self.prefix == "":
            self.prefix = "phonons"

        if len(fixatoms_dof) > 0:
            self.fixatoms_dof = fixatoms_dof
            if self.mode == "enmfd" or self.mode == "nmfd":
                raise ValueError("Fixatoms is not implemented for the selected mode.")
        else:
            self.fixatoms_dof = np.array([])

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        super(DynMatrixMover, self).bind(ens, beads, nm, cell, bforce, prng, omaker)

        # Raises error for nbeads not equal to 1.
        if self.beads.nbeads > 1:
            raise ValueError(
                "Calculation not possible for number of beads greater than one."
            )

        self.ism = 1 / xp.sqrt(self.beads.m3[-1])
        self.phcalc.bind(self)

        self.dbeads = self.beads.clone()
        self.dcell = self.cell.clone()
        self.dforces = self.forces.clone(self.dbeads, self.dcell)

    def step(self, step=None):
        """Executes one step of phonon computation."""
        if step < 3 * self.beads.natoms:
            self.phcalc.step(step)
        else:
            self.phcalc.transform()
            self.refdynmatrix = apply_asr(self.asr, self.refdynmatrix, self.beads)
            self.printall(
                self.prefix, self.refdynmatrix, fixatoms_dof=self.fixatoms_dof
            )
            self.finish(
                status="success",
                message="Dynamic matrix is calculated. Exiting simulation",
            )

    def printall(self, prefix, dmatx, deltaw=0.0, fixatoms_dof=np.array([])):
        """Prints out diagnostics for a given dynamical matrix."""

        dmatx = dmatx + xp.eye(len(dmatx)) * deltaw
        if deltaw != 0.0:
            wstr = " !! Shifted by %e !!" % (deltaw)
        else:
            wstr = ""

        # Get active arrays:
        activedof = 3 * self.beads.natoms - fixatoms_dof.size
        if fixatoms_dof.size > 0:
            active_atoms_mask = np.delete(
                np.arange(3 * self.beads.natoms), fixatoms_dof
            )
        else:
            active_atoms_mask = np.arange(3 * self.beads.natoms)

        dmatx_full = dmatx
        ism_full = self.ism
        dmatx = dmatx[active_atoms_mask][:, active_atoms_mask]
        ism = self.ism[active_atoms_mask]

        # prints out the dynamical matrix
        outfile = self.output_maker.get_output(self.prefix + ".dynmat", "w")
        outfile.write("# Dynamical matrix (atomic units)" + wstr + "\n")
        for i in range(activedof):
            outfile.write(" ".join(map(str, to_numpy(dmatx[i]))) + "\n")
        outfile.close_stream()

        # prints out the Hessian for the activedof
        outfile = self.output_maker.get_output(self.prefix + ".hess", "w")
        outfile.write("# Hessian matrix (atomic units)" + wstr + "\n")
        for i in range(activedof):
            outfile.write(
                " ".join(map(str, to_numpy(dmatx[i] / (ism[i] * ism)))) + "\n"
            )
        outfile.close_stream()

        # prints out the full Hessian (with all zeros)
        outfile = self.output_maker.get_output(self.prefix + "_full.hess", "w")
        outfile.write("# Hessian matrix (atomic units)" + wstr + "\n")
        for i in range(3 * self.beads.natoms):
            outfile.write(
                " ".join(map(str, to_numpy(dmatx_full[i] / (ism_full[i] * ism_full))))
                + "\n"
            )
        outfile.close_stream()

        eigsys = xp.linalg.eigh(dmatx)

        # prints eigenvalues
        outfile = self.output_maker.get_output(self.prefix + ".eigval", "w")
        outfile.write("# Eigenvalues (atomic units)" + wstr + "\n")
        outfile.write("\n".join(map(str, to_numpy(eigsys[0]))))
        outfile.close_stream()

        # prints eigenvectors
        outfile = self.output_maker.get_output(self.prefix + ".eigvec", "w")
        outfile.write(
            "# Eigenvector  matrix from the dynamical matrix (normalized)" + "\n"
        )
        for i in range(activedof):
            outfile.write(" ".join(map(str, to_numpy(eigsys[1][i]))) + "\n")
        outfile.close_stream()

        # prints eigenmodes
        eigmode = 1.0 * eigsys[1]
        for i in range(activedof):
            eigmode[i] *= ism[i]
        for i in range(activedof):
            eigmode[:, i] /= xp.sqrt(eigmode[:, i] @ eigmode[:, i])
        outfile = self.output_maker.get_output(self.prefix + ".mode", "w")

        outfile.write("# Phonon modes (cartesian space and normalized)" + "\n")
        for i in range(activedof):
            outfile.write(" ".join(map(str, to_numpy(eigmode[i]))) + "\n")
        outfile.close_stream()


class DummyPhononCalculator:
    """No-op PhononCalculator"""

    def __init__(self):
        pass

    def bind(self, dm):
        """Reference all the variables for simpler access."""
        self.dm = dm

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""
        pass

    def transform(self):
        """Dummy transformation step which does nothing."""
        pass


class FDPhononCalculator(DummyPhononCalculator):
    """Finite difference phonon evaluator."""

    def bind(self, dm):
        """Reference all the variables for simpler access."""
        super(FDPhononCalculator, self).bind(dm)

        nq = xp_size(self.dm.beads.q)

        # the dynamical matrix and the refined one, of size 3 natoms x 3 natoms
        self.dm.dynmatrix = xp.asarray(
            zeros_or_check(self.dm.dynmatrix, (nq, nq), "Force constant matrix")
        ).reshape((nq, nq))
        self.dm.refdynmatrix = xp.asarray(
            zeros_or_check(
                self.dm.refdynmatrix, (nq, nq), "Refined force constant matrix"
            )
        ).reshape((nq, nq))

    def step(self, step=None):
        """Computes one row of the dynamic matrix."""

        if step not in self.dm.fixatoms_dof:
            # initializes the finite deviation
            dev = xp.zeros(3 * self.dm.beads.natoms)
            dev[step] = self.dm.deltax
            # displaces kth d.o.f by delta.
            self.dm.dbeads.q = self.dm.beads.q + dev
            plus = -self.dm.dforces.f
            # displaces kth d.o.f by -delta.
            self.dm.dbeads.q = self.dm.beads.q - dev
            minus = -self.dm.dforces.f
            # computes a row of force-constant matrix
            dmrow = (
                (plus - minus) / (2 * self.dm.deltax) * self.dm.ism[step] * self.dm.ism
            )
            self.dm.dynmatrix[step] = dmrow
            self.dm.refdynmatrix[step] = dmrow
        else:
            info(" We have skipped the dof # {}.".format(step), verbosity.low)

    def transform(self):
        dm = self.dm.dynmatrix
        rdm = self.dm.dynmatrix
        self.dm.dynmatrix = 0.50 * (dm + dm.T)
        self.dm.refdynmatrix = 0.50 * (rdm + rdm.T)


class NMFDPhononCalculator(FDPhononCalculator):
    """Normal mode finite difference phonon evaluator."""

    def bind(self, dm):
        """Reference all the variables for simpler access."""
        super(NMFDPhononCalculator, self).bind(dm)

        if bool(xp.all(self.dm.dynmatrix == 0)):
            raise ValueError("Force constant matrix size not found")

        # the refined dynamical matrix, of size 3 natoms x 3 natoms
        nq = xp_size(self.dm.beads.q)
        self.dm.refdynmatrix = zeros_or_check(
            self.dm.refdynmatrix, (nq, nq), "Refined force constant matrix"
        )

        self.dm.w2, self.dm.U = xp.linalg.eigh(self.dm.dynmatrix)
        # scales each row by ism
        self.dm.V = self.dm.U * self.dm.ism[:, None]

    def step(self, step=None):
        """Computes one row of the dynamic matrix."""

        # initializes the finite deviation
        vknorm = xp.sqrt(self.dm.V[:, step] @ self.dm.V[:, step])
        dev = xp.real(self.dm.V[:, step] / vknorm) * self.dm.deltax
        # displaces by -delta along kth normal mode.
        self.dm.dbeads.q = self.dm.beads.q + dev
        plus = -self.dm.dforces.f.reshape(-1)
        # displaces by -delta along kth normal mode.
        self.dm.dbeads.q = self.dm.beads.q - dev
        minus = -self.dm.dforces.f.reshape(-1)
        # computes a row of the refined dynmatrix, in the basis of the eigenvectors of the first dynmatrix
        dmrowk = (plus - minus) / (2 * self.dm.deltax / vknorm)
        self.dm.refdynmatrix[step] = self.dm.V.T @ dmrowk

    def transform(self):
        self.dm.refdynmatrix = self.dm.U @ (self.dm.refdynmatrix @ self.dm.U.T)
        rdm = self.dm.refdynmatrix
        self.dm.refdynmatrix = 0.50 * (rdm + rdm.T)


class ENMFDPhononCalculator(NMFDPhononCalculator):
    """Energy scaled normal mode finite difference phonon evaluator."""

    def step(self, step=None):
        """Computes one row of the dynamic matrix."""

        # initializes the finite deviation
        vknorm = xp.sqrt(self.dm.V[:, step] @ self.dm.V[:, step])
        edelta = vknorm * xp.sqrt(self.dm.deltae * 2.0 / abs(self.dm.w2[step]))
        if edelta > 100 * self.dm.deltax:
            edelta = 100 * self.dm.deltax
        dev = xp.real(self.dm.V[:, step] / vknorm) * edelta
        # displaces by -delta along kth normal mode.
        self.dm.dbeads.q = self.dm.beads.q + dev
        plus = -self.dm.dforces.f.reshape(-1)
        # displaces by -delta along kth normal mode.
        self.dm.dbeads.q = self.dm.beads.q - dev
        minus = -self.dm.dforces.f.reshape(-1)
        # computes a row of the refined dynmatrix, in the basis of the eigenvectors of the first dynmatrix
        dmrowk = (plus - minus) / (2 * edelta / vknorm)
        self.dm.refdynmatrix[step] = self.dm.V.T @ dmrowk
