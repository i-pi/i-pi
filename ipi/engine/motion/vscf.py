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

__all__ = ["NormalModeMover"]

import numpy as np
import os

from ipi.engine.motion import Motion
from ipi.utils.depend import *

# from ipi.utils import units
from ipi.utils.phonontools import apply_asr
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, info

# from ipi.utils.io import print_file
# from ipi.engine.atoms import Atoms

try:
    import scipy
    from scipy.special import hermite
    from scipy.special import logsumexp
    from scipy.interpolate import interp1d
    from scipy.interpolate import interp2d
except Exception as e:
    scipy = None
    scipy_exception = e
from itertools import combinations


class NormalModeMover(Motion):
    """Normal Mode analysis."""

    def __init__(
        self,
        fixcom=False,
        fixatoms=None,
        mode="imf",
        dynmat=np.zeros(0, float),
        prefix="",
        asr="none",
        nprim="1",
        fnmrms="1.0",
        nevib="25.0",
        nint="101",
        nbasis="10",
        athresh="1e-2",
        ethresh="1e-2",
        nkbt="4.0",
        nexc="5",
        solve=False,
        grid=True,
        mptwo=False,
        print_mftpot=False,
        print_1b_map=False,
        print_2b_map=False,
        threebody=False,
        print_vib_density=False,
        nparallel=1,
        alpha=1.0,
        pair_range=np.zeros(0, int),
    ):
        """Initialises NormalModeMover.
        Args:
        fixcom : An optional boolean which decides whether the centre of mass
                  motion will be constrained or not. Defaults to False.
        dynmatrix : A 3Nx3N array that stores the dynamic matrix.
        refdynmatrix : A 3Nx3N array that stores the refined dynamic matrix.
        """

        super(NormalModeMover, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        # Finite difference option.
        self.mode = mode
        if self.mode == "imf":
            self.calc = IMF()
        elif self.mode == "vscf":
            self.calc = VSCF()

        self.dynmatrix = dynmat
        self.mode = mode
        self.frefine = False
        self.U = None
        self.V = None
        self.prefix = prefix
        self.asr = asr
        self.nprim = nprim  # 1
        self.nz = 0  # 1
        self.alpha = alpha  # 1
        self.fnmrms = fnmrms  # 1.0
        self.nevib = nevib  # 25.0
        self.nint = nint  # 101
        self.nbasis = nbasis  # 10
        self.athresh = athresh  # 1e-2
        self.ethresh = ethresh  # 1e-2
        self.nkbt = nkbt  # 4.0
        self.nexc = nexc  # 5
        self.solve = solve  # 5
        self.grid = grid  # 5
        self.mptwo = mptwo  # False
        self.print_mftpot = print_mftpot  # False
        self.print_1b_map = print_1b_map  # False
        self.print_2b_map = print_2b_map  # False
        self.threebody = threebody  # False
        self.print_vib_density = print_vib_density  # False
        self.nparallel = nparallel  # False
        self.pair_range = pair_range

        if self.prefix == "":
            self.prefix = self.mode

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        super(NormalModeMover, self).bind(ens, beads, nm, cell, bforce, prng, omaker)

        # TODO - A LOT OF THESE DEFINITIONS ARE NOT USING DEPEND OBJECTS CORRECTLY
        # MANY OF THESE SHOUD BE LINKED THROUGH DPIPE, NOT DEREFERENCED. THIS WOULD
        # BREAK IF USED WITH REPLICA EXCHANGE FOR INSTANCE.

        self.temp = self.ensemble.temp

        # Raises error for nbeads not equal to 1.
        if self.beads.nbeads > 1:
            raise ValueError(
                "Calculation not possible for number of beads greater than one."
            )

        self.ism = 1 / np.sqrt(dstrip(self.beads.m3[-1]))
        self.m = dstrip(self.beads.m)
        self.calc.bind(self)

        self.dbeads = self.beads.copy(nbeads=self.nparallel)
        self.dcell = self.cell.copy()
        self.dforces = self.forces.copy(self.dbeads, self.dcell)

    def step(self, step=None):
        """Executes one step of phonon computation."""
        self.calc.step(step)


class DummyCalculator:
    """No-op Calculator"""

    def __init__(self):
        pass

    def bind(self, imm):
        """Reference all the variables for simpler access."""
        self.imm = imm

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""
        pass


class IMF(DummyCalculator):
    """Temperature scaled normal mode Born-Oppenheimer surface evaluator."""

    def bind(self, imm):
        """Reference all the variables for simpler access."""

        if scipy is None:
            info(" @NM: scipy import failed", verbosity.low)
            raise scipy_exception

        super(IMF, self).bind(imm)

        # Initialises a 3*number of atoms X 3*number of atoms dynamic matrix.
        if self.imm.dynmatrix.size != (self.imm.beads.q.size * self.imm.beads.q.size):
            if self.imm.dynmatrix.size == 0:
                self.imm.dynmatrix = np.zeros(
                    (self.imm.beads.q.size, self.imm.beads.q.size), float
                )
            else:
                raise ValueError(
                    "Force constant matrix size does not match system size"
                )
        else:
            self.imm.dynmatrix = self.imm.dynmatrix.reshape(
                ((self.imm.beads.q.size, self.imm.beads.q.size))
            )

        # Applies ASR
        self.imm.dynmatrix = apply_asr(self.imm.asr, self.imm.dynmatrix, self.imm.beads)

        # Calculates normal mode directions.
        self.imm.w2, self.imm.U = np.linalg.eigh(self.imm.dynmatrix)
        # Sets
        self.imm.w2[np.abs(self.imm.w2) < 1e-15] = 0.0
        self.imm.nz = sum(np.abs(self.imm.w2) < 1e-15)

        # Calculates the normal mode frequencies.
        # TODO : Should also handle soft modes.
        self.imm.w = self.imm.w2 * 0
        self.imm.w[self.imm.nz :] = np.sqrt(self.imm.w2[self.imm.nz :])

        self.imm.V = self.imm.U.copy()
        for i in range(len(self.imm.V)):
            self.imm.V[:, i] *= self.imm.ism

        # Harm ZP RMS displacement along normal mode
        # Not temperature dependent so that sampled potentials can easily be
        # reused to evaluate free energy at different temperature.
        # self.imm.nmrms = np.zeros(len(self.imm.w))
        # self.imm.nmrms[self.imm.nz:] = np.sqrt( 0.5 / self.imm.w[self.imm.nz:])

        # Uses a temp dependent thermal RMS displacement as the scale of displacements.
        # TODO : Add a separate mapping temperature.
        self.imm.nmrms = np.zeros(len(self.imm.w))
        self.imm.nmrms[self.imm.nz :] = np.sqrt(
            0.5
            / self.imm.w[self.imm.nz :]
            / np.tanh(self.imm.w[self.imm.nz :] / self.imm.temp / 2.0)
        )
        self.nmrms = self.imm.nmrms

        # Similarly this is also not temperature dependent.
        # self.imm.nmevib =  0.5 * self.imm.w

        # Harm vibr energy at finite temp
        self.imm.nmevib = np.zeros(len(self.imm.w))
        self.imm.nmevib[self.imm.nz :] = (
            0.50
            * self.imm.w[self.imm.nz :]
            / np.tanh(self.imm.w[self.imm.nz :] / self.imm.temp / 2.0)
        )

        # Fraction of the harmonic RMS displacement
        # used to sample along a normal mode
        self.fnmrms = self.imm.fnmrms

        # Multiple of harmonic vibrational energy up to which
        # the BO surface is sampled
        self.nevib = self.imm.nevib

        # Integration points for numerical integration of
        # Hamiltonian matrix elements
        self.nint = self.imm.nint

        # Number of SHO states used as basis for anharmonic wvfn
        self.nbasis = self.imm.nbasis

        # Convergence threshold for fractional error in vibrational free energy
        self.athresh = self.imm.athresh

        # Number of primitive unit (cells) per simulation (cell)
        self.nprim = self.imm.nprim

        # Total number of vibrational modes
        self.dof = 3 * self.imm.beads.natoms

        # Initializes the (an)har energy correction
        self.total_anhar_free_energy = 0.0
        self.total_har_free_energy = 0.0
        self.total_anhar_internal_energy = 0.0
        self.total_har_internal_energy = 0.0

        # Potential energy (offset) at equilibrium positions per primitve unit (cell)
        self.v0 = 0

        # Initializes the SHO wvfn basis for solving the 1D Schroedinger equation
        self.hermite_functions = [hermite(i) for i in range(max(20, 4 * self.nbasis))]

        # Sets the total number of steps for IMF.
        self.total_steps = 3 * self.imm.beads.natoms

    def psi(self, n, m, hw, q):
        """
        Returns the value of the n^th wavefunction of a
        harmonic oscillator with mass m, frequency freq
        at position x.
        """

        # Defines variables for easier referencing
        alpha = m * hw
        try:
            herfun = self.hermite_functions[n]
        except:
            herfun = hermite(n)

        norm = (alpha / np.pi) ** 0.25 * np.sqrt(
            1.0 / (2.0 ** (n) * np.math.factorial(n))
        )
        psival = norm * np.exp(-alpha * q**2 / 2.0) * herfun(np.sqrt(alpha) * q)

        return psival

    def solve_schroedingers_equation(self, hw, psigrid, vgrid, return_eigsys=False):
        """
        Given the frequency of the HO, the wavefunctions and
        the gridded potential, solves the schrodinger's
        equation.
        """

        # The number of basis set elements.
        nbasis = len(psigrid)

        # Constructs the Hamiltonian matrix.
        h = np.zeros((nbasis, nbasis))
        for i in range(nbasis):
            for j in range(i, nbasis, 1):
                h[i][j] = np.dot(psigrid[i] * psigrid[j], vgrid)
            h[i][i] *= 0.5
            h[i][i] += 0.5 * (i + 0.5) * hw
        h += h.T

        # Diagonalise Hamiltonian matrix and evaluate anharmonic free energy and vibrational freq
        evals, evecs = np.linalg.eigh(h)

        # Calculates the free and internal energy
        A = -logsumexp(-1.0 * evals / self.imm.temp) * self.imm.temp
        E = np.sum(evals * np.exp(-1.0 * evals / self.imm.temp)) / np.sum(
            np.exp(-1.0 * evals / self.imm.temp)
        )

        if return_eigsys:
            return A, E, evals, evecs.T
        else:
            return A, E

    def step(self, step=None):
        """Computes the Born Oppenheimer curve along a normal mode."""

        if step == self.total_steps:
            self.terminate()

        # Ignores (near) zero modes.
        if step < self.imm.nz:
            info(" @NM: Ignoring the zero mode.", verbosity.medium)
            info(" ", verbosity.medium)
            return
        elif self.imm.w[step] < 9.1126705e-06:
            info(
                " @NM: Ignoring normal mode no.  %8d with frequency %15.8f cm^-1."
                % (step, self.imm.w[step] * 219474),
                verbosity.medium,
            )
            info(" ", verbosity.medium)
            self.imm.nz += 1
            return
        else:
            info(
                " @NM: Treating normal mode no.  %8d with frequency %15.8f cm^-1."
                % (step, self.imm.w[step] * 219474),
                verbosity.medium,
            )

        self.v_indep_filename = (
            self.imm.output_maker.prefix
            + "."
            + self.imm.prefix
            + "."
            + str(step)
            + ".qvf"
        )
        if os.path.exists(self.v_indep_filename):
            # Loads the displacemnts, the potential energy and the forces.
            qlist, vlist, flist = np.loadtxt(self.v_indep_filename).T

            # Fits cubic splines to data.
            info(" @NM: Fitting cubic splines.", verbosity.medium)
            vspline = interp1d(qlist, vlist, kind="cubic", bounds_error=False)

            # Converge wrt size of SHO basis
            bs_iter = 0
            bs_Aanh = [1e-20]
            bs_Eanh = [1e-20]

            # Calculates the 1D correction potential on a grid.
            qgrid = np.linspace(np.min(qlist), np.max(qlist), self.nint)
            vgrid = np.asarray([(np.nan_to_num(vspline(x))) for x in qgrid])

            while True:
                nnbasis = max(1, self.nbasis - 5) + 5 * bs_iter
                nnbasis = self.nbasis + 5 * bs_iter

                psigrid = np.zeros((nnbasis, self.nint))
                # Calculates the wavefunctions on the grid.
                for i in range(nnbasis):
                    psigrid[i] = self.psi(i, 1.0, self.imm.w[step], qgrid)
                    psigrid[i] = psigrid[i] / np.sqrt(np.dot(psigrid[i], psigrid[i]))

                # Solves the Schroedinger's Equation.
                bs_AEanh = self.solve_schroedingers_equation(
                    self.imm.w[step], psigrid, vgrid
                )
                bs_Aanh.append(bs_AEanh[0])
                bs_Eanh.append(bs_AEanh[1])

                dA = np.abs(bs_Aanh[-1] - bs_Aanh[-2]) / (self.dof - self.imm.nz)
                info(
                    " @NM: CONVERGENCE : nbasis = %5d    A =  %10.8e   D(A) =  %10.8e /  %10.8e"
                    % (nnbasis, bs_Aanh[-1], dA, self.athresh),
                    verbosity.medium,
                )

                # Check whether anharmonic frequency is converged
                if dA < self.athresh:
                    break

                bs_iter += 1
            Aanh = [bs_Aanh[-1]]
            Eanh = [bs_Eanh[-1]]

        else:
            # initializes the list containing sampled displacements,
            # forces and energies.
            vlist = []
            flist = []
            qlist = []

            # Adds to the list of "sampled configuration" the one that
            # corresponds to the minimum.
            q0 = 0.0
            v0 = dstrip(self.imm.forces.pots).copy()[0] / self.nprim
            f0 = (
                np.dot(dstrip(self.imm.forces.f).copy()[0], np.real(self.imm.V.T[step]))
                / self.nprim
            )

            self.v0 = v0  # TODO CHECK IF NECESSARY

            vlist.append(0.0)
            flist.append(f0)
            qlist.append(q0)

            # Converge anharmonic vibrational energies w.r.t. density of sampling points
            Aanh = []
            Eanh = []
            Aanh.append(1e-20)

            sampling_density_iter = -1

            while True:
                sampling_density_iter += 1

                # Doubles the grid spacing, so that an estimate of the
                # anharmonic free energy convergence is
                # possible at default/input grid spacing
                ffnmrms = self.fnmrms * 0.5**sampling_density_iter * 2.0

                # Calculates the displacement in Cartesian coordinates.
                nmd = ffnmrms * self.imm.nmrms[step]
                dev = np.real(self.imm.V.T[step]) * nmd * np.sqrt(self.nprim)

                # After the first iteration doubles the displacement to avoid
                # calculation of the potential at configurations already v_indep_listited
                # in the previous iteration.
                if sampling_density_iter == 0:
                    delta_counter = 1
                else:
                    delta_counter = 2

                counter = 1

                # Explores configurations until the sampled energy exceeds
                # a user-defined threshold of the zero-point energy.
                while True:
                    # Displaces along the normal mode.
                    self.imm.dbeads.q = self.imm.beads.q + dev * counter

                    # Stores the "anharmonic" component of the potential
                    # and the force.
                    dv = (
                        dstrip(self.imm.dforces.pots).copy()[0] / self.nprim
                        - 0.50 * self.imm.w2[step] * (nmd * counter) ** 2
                        - v0
                    )
                    df = np.dot(
                        dstrip(self.imm.dforces.f).copy()[0],
                        np.real(self.imm.V.T[step]),
                    ) / self.nprim + self.imm.w2[step] * (nmd * counter)

                    # Adds to the list.
                    # Also stores the total energetics i.e. including
                    # the harmonic component.
                    vlist.append(dv)
                    flist.append(df)
                    qlist.append(nmd * counter)

                    # Bailout condition.
                    if self.nevib * self.imm.nmevib[step] < np.abs(
                        0.50 * self.imm.w2[step] * (nmd * counter) ** 2 + dv
                    ):
                        break

                    # Increases the displacement by 1 or 2 depending on the iteration.
                    counter += delta_counter

                info(
                    " @NM: Using %8d configurations along the +ve direction."
                    % (counter,),
                    verbosity.medium,
                )

                counter = -1

                # Similarly displaces in the "-ve direction"
                while True:
                    # Displaces along the normal mode.
                    self.imm.dbeads.q = self.imm.beads.q + dev * counter

                    # Stores the "anharmonic" component of the potential
                    # and the force.
                    dv = (
                        dstrip(self.imm.dforces.pots).copy()[0] / self.nprim
                        - 0.50 * self.imm.w2[step] * (nmd * counter) ** 2
                        - v0
                    )
                    df = np.dot(
                        dstrip(self.imm.dforces.f).copy()[0],
                        np.real(self.imm.V.T[step]),
                    ) / self.nprim + self.imm.w2[step] * (nmd * counter)

                    # Adds to the list.
                    # Also stores the total energetics i.e. including
                    # the harmonic component.
                    vlist.append(dv)
                    flist.append(df)
                    qlist.append(nmd * counter)

                    # Bailout condition.
                    if self.nevib * self.imm.nmevib[step] < np.abs(
                        0.50 * self.imm.w2[step] * (nmd * counter) ** 2 + dv
                    ):
                        break

                    # Increases the displacement by 1 or 2 depending on the iteration.
                    counter -= delta_counter

                info(
                    " @NM: Using %8d configurations along the -ve direction."
                    % (-counter,),
                    verbosity.medium,
                )

                # Fits cubic splines to data.
                info("@NM: Fitting cubic splines.", verbosity.medium)
                vspline = interp1d(qlist, vlist, kind="cubic", bounds_error=False)

                # Converge wrt size of SHO basis
                bs_iter = 0
                bs_Aanh = [1e-20]
                bs_Eanh = [1e-20]

                # Calculates the 1D correction potential on a grid.
                qgrid = np.linspace(np.min(qlist), np.max(qlist), self.nint)
                vgrid = np.asarray([(np.nan_to_num(vspline(x))) for x in qgrid])

                while True:
                    nnbasis = max(1, self.nbasis - 5) + 5 * bs_iter

                    psigrid = np.zeros((nnbasis, self.nint))
                    # Calculates the wavefunctions on the grid.
                    for i in range(nnbasis):
                        psigrid[i] = self.psi(i, 1.0, self.imm.w[step], qgrid)
                        psigrid[i] = psigrid[i] / np.sqrt(
                            np.dot(psigrid[i], psigrid[i])
                        )

                    # Solves the Schroedinger's Equation.
                    bs_AEanh = self.solve_schroedingers_equation(
                        self.imm.w[step], psigrid, vgrid
                    )
                    bs_Aanh.append(bs_AEanh[0])
                    bs_Eanh.append(bs_AEanh[1])

                    dA = np.abs(bs_Aanh[-1] - bs_Aanh[-2]) / (self.dof - self.imm.nz)
                    info(
                        " @NM: CONVERGENCE : fnmrms = %10.8e   nbasis = %5d    A =  %10.8e   D(A) =  %10.8e /  %10.8e"
                        % (ffnmrms, nnbasis, bs_Aanh[-1], dA, self.athresh),
                        verbosity.medium,
                    )

                    # Check whether anharmonic frequency is converged
                    if dA < self.athresh:
                        break

                    bs_iter += 1

                Aanh.append(bs_Aanh[-1])
                Eanh.append(bs_Eanh[-1])

                # Check whether anharmonic frequency is converged
                dA = np.abs(Aanh[-1] - Aanh[-2]) / (self.dof - self.imm.nz)
                if dA < self.athresh:
                    break

            # Prints the normal mode displacement, the potential and the force.
            outfile = self.imm.output_maker.get_output(
                self.imm.prefix + "." + str(step) + ".qvf"
            )
            np.savetxt(
                outfile,
                np.c_[qlist, vlist, flist],
                header="Frequency = %10.8f" % self.imm.w[step],
            )
            outfile.close_stream()

        # prints the mapped potential.
        if self.imm.print_1b_map is True:
            output_grid = np.linspace(np.min(qlist), np.max(qlist), 100)
            outfile = self.imm.output_maker.get_output(
                self.imm.prefix + "." + str(step) + ".vfit"
            )
            np.savetxt(
                outfile,
                np.c_[output_grid, vspline(output_grid)],
                header="Frequency = %10.8f" % self.imm.w[step],
            )
            info(
                " @NM: Prints the mapped potential energy to %s"
                % (self.imm.prefix + "." + str(step) + ".vfit"),
                verbosity.medium,
            )
            outfile.close_stream()

        # Done converging wrt size of SHO basis.
        # Calculates the harmonic free and internal energy.
        Ahar = (
            -logsumexp(
                [
                    -1.0 * np.sqrt(self.imm.w2[step]) * (0.5 + i) / self.imm.temp
                    for i in range(nnbasis)
                ]
            )
            * self.imm.temp
        )
        Zhar = np.sum(
            [
                np.exp(-1.0 * np.sqrt(self.imm.w2[step]) * (0.5 + i) / self.imm.temp)
                for i in range(nnbasis)
            ]
        )
        Ehar = (
            np.sum(
                [
                    np.sqrt(self.imm.w2[step])
                    * (0.5 + i)
                    * np.exp(
                        -1.0 * np.sqrt(self.imm.w2[step]) * (0.5 + i) / self.imm.temp
                    )
                    for i in range(nnbasis)
                ]
            )
            / Zhar
        )

        info(
            " @NM: HAR frequency     =  %10.8e" % (self.imm.w[step],), verbosity.medium
        )
        info(" @NM: HAR free energy   =  %10.8e" % (Ahar,), verbosity.medium)
        info(" @NM: IMF free energy   =  %10.8e" % (Aanh[-1],), verbosity.medium)
        info("\n", verbosity.medium)
        self.total_anhar_free_energy += Aanh[-1]
        self.total_har_free_energy += Ahar
        self.total_anhar_internal_energy += Eanh[-1]
        self.total_har_internal_energy += Ehar

    def terminate(self):
        """
        Prints out the free and internal energy
        for HAR and IMF, and triggers a soft exit.
        """

        info(
            " @NM: Potential offset               =  %10.8e" % (self.v0,), verbosity.low
        )
        info(
            " @NM: HAR free energy                =  %10.8e"
            % (
                np.sum(
                    (
                        0.5 * np.sqrt(self.imm.w2[self.imm.nz :])
                        + self.imm.temp
                        * np.log(
                            1.0
                            - np.exp(
                                -np.sqrt(self.imm.w2[self.imm.nz :]) / self.imm.temp
                            )
                        )
                    )
                )
                / self.nprim
                + self.v0,
            ),
            verbosity.low,
        )
        info(
            " @NM: IMF free energy correction     =  %10.8e"
            % (
                (self.total_anhar_free_energy - self.total_har_free_energy)
                / self.nprim,
            ),
            verbosity.low,
        )
        info(
            " @NM: HAR internal energy            =  %10.8e"
            % (
                np.sum(
                    np.sqrt(self.imm.w2[self.imm.nz :])
                    * (
                        0.5
                        + 1.0
                        / (
                            np.exp(np.sqrt(self.imm.w2[self.imm.nz :]) / self.imm.temp)
                            - 1
                        )
                    )
                )
                / self.nprim
                + self.v0,
            ),
            verbosity.low,
        )
        info(
            " @NM: IMF internal energy correction =  %10.8e"
            % (
                (self.total_anhar_internal_energy - self.total_har_internal_energy)
                / self.nprim,
            ),
            verbosity.low,
        )
        info(
            " @NM: ALL QUANTITIES PER PRIMITIVE UNIT CELL (WHERE APPLICABLE) \n",
            verbosity.low,
        )
        softexit.trigger(
            status="success", message=" @NM: The IMF calculation has terminated."
        )


class VSCF(IMF):
    """ """

    def bind(self, imm):
        """
        Reference all the variables for simpler access
        """

        if scipy is None:
            info(" @NM: scipy import failed", verbosity.low)
            raise scipy_exception

        super(VSCF, self).bind(imm)
        self.nz = self.imm.nz
        # self.print_2b_map = self.imm.print_2b_map
        # self.threebody = self.imm.threebody
        self.solve = self.imm.solve
        self.grid = self.imm.grid
        self.alpha = self.imm.alpha
        self.pair_range = self.imm.pair_range

        # Filenames for storing the number of samples configurations
        # and their potential energy.
        self.modes_filename = "modes.dat"
        self.v_offset_prefix = "potoffset"
        self.npts_pos_prefix = "npts_pos"
        self.npts_neg_prefix = "npts_neg"
        self.v_indep_file_prefix = "vindep"
        self.v_indep_grid_file_prefix = "vindep_grid"
        self.v_coupled_file_prefix = "vcoupled"
        self.v_coupled_grid_file_prefix = "vcoupled_grid"
        self.v_coupled_grids_file_prefix = "vcoupled_grids"
        self.eigenvalue_filename = "evals.dat"
        self.eigenvector_filename = "evecs.dat"

        # Creates a list of modes with frequencies greater than 2 cm^-1.
        if os.path.exists(self.imm.output_maker.prefix + "." + self.modes_filename):
            self.inms = np.loadtxt(
                self.imm.output_maker.prefix + "." + self.modes_filename, dtype=int
            ).tolist()
        else:
            info(" @NM: Identifying relevant frequency modes.", verbosity.medium)
            self.inms = []
            for inm in range(self.dof):
                if self.imm.w[inm] < 9.1126705e-06:
                    info(
                        " @NM: Ignoring normal mode no.  %8d with frequency %15.8f cm^-1."
                        % (
                            inm,
                            self.imm.w[inm] * 219474,
                        ),
                        verbosity.medium,
                    )
                    continue
                else:
                    self.inms.append(inm)

            # Save for use in VSCFSolver.
            outfile = self.imm.output_maker.get_output(self.modes_filename)
            np.savetxt(
                outfile,
                self.inms,
                fmt="%i",
                header="Indices of modes that are considered in the calculation.",
            )
            outfile.close_stream()

        # Saves the total number of steps for automatic termination.
        dof = len(self.inms)
        self.total_steps = dof + dof * (dof - 1) / 2

        # Saves the indices of pairs of modes
        self.pair_combinations = list(combinations(self.inms, 2))

        # Selects the range of pair of modes to be calculated.
        if self.pair_range.size == 0:
            self.pair_range = np.asarray([0, len(self.pair_combinations)])

        # Variables for storing the number of sampled configurations
        # along the +ve and -ve displacements along normal modes and
        # the sampled potential energy.
        self.npts = np.zeros(self.dof, dtype=int)
        self.npts_neg = np.zeros(self.dof, dtype=int)
        self.npts_pos = np.zeros(self.dof, dtype=int)
        self.v_indep_list = []
        self.displacements_nm = []
        self.displacements_nm = []

        for x in range(self.dof - dof):
            self.displacements_nm.append([0])

        if self.solve or self.grid:
            self.q_grids = np.zeros((self.dof, self.nint))
            self.v_indep_grids = np.zeros((self.dof, self.nint))
            self.v_mft_grids = np.zeros((self.dof, self.nint))

            self.psi_i_grids = np.zeros((self.dof, self.nbasis, self.nint))
            self.rho_grids = np.zeros((self.dof, self.nint))

            self.evecs_imf = np.zeros((self.dof, self.nbasis, self.nbasis))
            self.evals_imf = np.zeros((self.dof, self.nbasis))

            self.evals_vscf = np.zeros((self.dof, self.nbasis))
            self.evecs_vscf = np.zeros((self.dof, self.nbasis, self.nbasis))

    def step(self, step=None):
        """Computes the Born-Oppenheimer curve along a normal mode."""

        # Performs some basic initialization.
        if step == 0:
            # Initialize overall potential offset
            self.v_offset_filename = self.v_offset_prefix + ".dat"
            if os.path.exists(
                self.imm.output_maker.prefix + "." + self.v_offset_filename
            ):
                self.v0 = np.loadtxt(
                    self.imm.output_maker.prefix + "." + self.v_offset_filename
                )
            else:
                self.v0 = dstrip(self.imm.forces.pots).copy()[0] / self.nprim
                outfile = self.imm.output_maker.get_output(self.v_offset_filename)
                np.savetxt(outfile, [self.v0])
                outfile.close_stream()

        # Maps 1D curves.
        elif step <= len(self.inms):
            # Selects the normal mode to map out.
            self.inm = self.inms[step - 1]

            # Defines the names of the files that store the sampled
            # and the interpolated potential energy.
            self.v_indep_filename = (
                self.v_indep_file_prefix + "." + str(self.inm) + ".dat"
            )
            self.v_indep_grid_filename = (
                self.v_indep_grid_file_prefix + "." + str(self.inm) + ".dat"
            )

            info(
                "\n @NM: Treating normal mode no.  %8d with frequency %15.8f cm^-1."
                % (self.inm, self.imm.w[self.inm] * 219474),
                verbosity.medium,
            )

            # If the indepent modes are already calculated, just loads from file.
            if os.path.exists(
                self.imm.output_maker.prefix + "." + self.v_indep_filename
            ):
                # Reads the number of sampled configurations from the header
                # and the sampeld potential energy from the body.
                with open(
                    self.imm.output_maker.prefix + "." + self.v_indep_filename
                ) as f:
                    header = [line.split() for line in f if line.startswith("#")][0]
                    self.npts_neg[self.inm] = int(header[2])
                    self.npts_pos[self.inm] = int(header[4])
                    self.npts[self.inms] = (
                        self.npts_neg[self.inms] + self.npts_pos[self.inms] + 1
                    )
                self.v_indep_list.append(
                    np.loadtxt(
                        self.imm.output_maker.prefix + "." + self.v_indep_filename
                    ).T
                )
                self.displacements_nm.append(
                    [
                        self.fnmrms
                        * self.nmrms[self.inm]
                        * (-self.npts_neg[self.inm] + i - 2.0)
                        for i in range(self.npts[self.inm] + 4)
                    ]
                )
                info(
                    " @NM: Loading the sampled potential energy for mode  %8d"
                    % (self.inm,),
                    verbosity.medium,
                )
                info(
                    " @NM: Using %8d configurations along the +ve direction."
                    % (self.npts_pos[self.inm],),
                    verbosity.medium,
                )
                info(
                    " @NM: Using %8d configurations along the -ve direction."
                    % (self.npts_neg[self.inm],),
                    verbosity.medium,
                )

            # If mapping has NOT been perfomed previously, maps the 1D curves.
            else:
                (
                    self.npts_neg[self.inm],
                    self.npts_pos[self.inm],
                    v_indeps,
                ) = self.one_dimensional_mapper(step)
                self.npts[self.inm] = (
                    self.npts_neg[self.inm] + self.npts_pos[self.inm] + 1
                )
                self.v_indep_list.append(v_indeps)
                self.displacements_nm.append(
                    [
                        self.fnmrms
                        * self.nmrms[self.inm]
                        * (-self.npts_neg[self.inm] + i - 2.0)
                        for i in range(self.npts[self.inm] + 4)
                    ]
                )

                info(
                    " @NM: Using %8d configurations along the +ve direction."
                    % (self.npts_pos[self.inm],),
                    verbosity.medium,
                )
                info(
                    " @NM: Using %8d configurations along the -ve direction."
                    % (self.npts_neg[self.inm],),
                    verbosity.medium,
                )
                info(
                    " @NM: Saving the sampled potential energy for mode  %8d in %s"
                    % (self.inm, self.v_indep_filename),
                    verbosity.medium,
                )
                outfile = self.imm.output_maker.get_output(self.v_indep_filename)
                np.savetxt(
                    outfile,
                    v_indeps,
                    header=" npts_neg: %10d npts_pos: %10d"
                    % (self.npts_neg[self.inm], self.npts_pos[self.inm]),
                )
                outfile.close_stream()

            # We need the independent mode correction on a grid.
            # Checks if the potential exists otherwise loads from file.
            if self.grid:
                if os.path.exists(
                    self.imm.output_maker.prefix + "." + self.v_indep_grid_filename
                ):
                    igrid, vigrid = np.loadtxt(
                        self.imm.output_maker.prefix + "." + self.v_indep_grid_filename
                    ).T
                else:
                    info(
                        " @NM: Interpolating the potential energy on a grid of %8d points."
                        % (self.nint,),
                        verbosity.medium,
                    )
                    vspline = interp1d(
                        self.displacements_nm[-1],
                        self.v_indep_list[-1],
                        kind="cubic",
                        bounds_error=False,
                    )
                    igrid = np.linspace(
                        -self.npts_neg[self.inm] * self.fnmrms * self.nmrms[self.inm],
                        self.npts_pos[self.inm] * self.fnmrms * self.nmrms[self.inm],
                        self.nint,
                    )
                    vigrid = np.asarray(
                        [
                            float(
                                vspline(igrid[iinm])
                                - 0.5 * self.imm.w2[self.inm] * igrid[iinm] ** 2
                                - self.v0
                            )
                            for iinm in range(self.nint)
                        ]
                    )

                    # Save coupling correction to file for vistualisation.
                    info(
                        " @NM: Saving the interpolated potential energy to %s"
                        % (self.v_indep_grid_filename,),
                        verbosity.medium,
                    )
                    outfile = self.imm.output_maker.get_output(
                        self.v_indep_grid_filename
                    )
                    np.savetxt(outfile, np.c_[igrid, vigrid])
                    outfile.close_stream()

                # Stores the interpolated potential in memory.
                self.q_grids[self.inm][:] = igrid
                self.v_indep_grids[self.inm][:] = vigrid

        # Maps 2D surfaces if the index of the pair lies in range.
        elif step - len(self.inms) - 1 < self.pair_range[1]:
            # Checks if the index lies in range.
            if not self.pair_range[0] <= step - len(self.inms) - 1:
                return

            # Selects the normal mode pair to map out.
            vijgrid = None
            self.inm, self.jnm = self.pair_combinations[step - len(self.inms) - 1]
            self.inm_index, self.jnm_index = self.inm - self.nz, self.jnm - self.nz

            # Defines the names of the files that store the sampled
            # and the interpolated potential energy.
            self.v_coupled_filename = (
                self.v_coupled_file_prefix
                + "."
                + str(self.inm)
                + "."
                + str(self.jnm)
                + ".dat"
            )
            self.v_coupled_grid_filename = (
                self.v_coupled_grid_file_prefix
                + "."
                + str(self.inm)
                + "."
                + str(self.jnm)
                + ".dat"
            )

            info(
                "\n @NM: Treating normal modes no.  %8d  and %8d  with frequencies %15.8f cm^-1 and %15.8f cm^-1, respectively."
                % (
                    self.inm,
                    self.jnm,
                    self.imm.w[self.inm] * 219474,
                    self.imm.w[self.jnm] * 219474,
                ),
                verbosity.medium,
            )

            # Skips the step if the grid file exists.
            if os.path.exists(
                self.imm.output_maker.prefix + "." + self.v_coupled_grid_filename
            ):
                return

            elif (
                os.path.exists(
                    self.imm.output_maker.prefix + "." + self.v_coupled_filename
                )
                is not True
            ):
                # Initializes the grid for interpolating the potential when
                # displacements are made along pairs of normal modes.
                self.v_coupled = np.zeros(
                    (self.npts[self.inm] + 4) * (self.npts[self.jnm] + 4)
                )

                # Calculates the displacements as linear combinations of displacements along independent modes.
                displacements_nmi = self.displacements_nm[self.inm]
                displacements_nmj = self.displacements_nm[self.jnm]

                # Calculates the potential energy at the displaced positions.
                k = 0
                didjv = []
                unit_displacement_nmi = np.real(self.imm.V.T[self.inm]) * np.sqrt(
                    self.nprim
                )
                unit_displacement_nmj = np.real(self.imm.V.T[self.jnm]) * np.sqrt(
                    self.nprim
                )
                info(
                    " @NM: Sampling a total of %8d configurations."
                    % (len(displacements_nmi) * len(displacements_nmj),),
                    verbosity.medium,
                )

                for i in range(self.npts[self.inm] + 4):
                    for j in range(self.npts[self.jnm] + 4):
                        # Uses on-axis potentials are available from 1D maps.
                        if (-self.npts_neg[self.inm] + i - 2) == 0:
                            self.v_coupled[k] = self.v_indep_list[self.jnm_index][j]
                        elif (-self.npts_neg[self.jnm] + j - 2) == 0:
                            self.v_coupled[k] = self.v_indep_list[self.inm_index][i]
                        else:
                            self.imm.dbeads.q = (
                                dstrip(self.imm.beads.q)
                                + displacements_nmi[i] * unit_displacement_nmi
                                + displacements_nmj[j] * unit_displacement_nmj
                            )
                            self.v_coupled[k] = (
                                dstrip(self.imm.dforces.pots)[0] / self.nprim
                            )
                        didjv.append(
                            [
                                displacements_nmi[i],
                                displacements_nmj[j],
                                self.v_coupled[k] - self.v0,
                            ]
                        )
                        k += 1

                # Saves the displacements and the sampled potential energy.
                info(
                    " @NM: Saving the sampled potential energy to %s."
                    % (self.v_coupled_filename,),
                    verbosity.medium,
                )
                outfile = self.imm.output_maker.get_output(self.v_coupled_filename)
                np.savetxt(outfile, didjv)
                outfile.close_stream()

            else:
                info(
                    " @NM: Skipping the mapping for modes %8d and %8d."
                    % (self.inm, self.jnm),
                    verbosity.medium,
                )
                if self.grid:
                    displacements_nmi, displacements_nmj, self.v_coupled = np.loadtxt(
                        self.imm.output_maker.prefix + "." + self.v_coupled_filename
                    ).T
                    self.v_coupled += self.v0
                    displacements_nmi = self.displacements_nm[self.inm]
                    displacements_nmj = self.displacements_nm[self.jnm]

            # We need the pair-wise coupling correction on a grid.
            # Checks if the correction exists otherwise loads from file.
            if self.grid:
                if os.path.exists(
                    self.imm.output_maker.prefix + "." + self.v_coupled_grid_filename
                ):
                    vijgrid = np.load(
                        self.imm.output_maker.prefix
                        + "."
                        + self.v_coupled_grid_filename
                    )

                else:
                    # Interpolates the displacements on a grid and saves for VSCFSOLVER.
                    info(
                        " @NM: Interpolating the potential energy on a %8d x %8d grid."
                        % (self.nint, self.nint),
                        verbosity.medium,
                    )
                    vtspl = interp2d(
                        displacements_nmi,
                        displacements_nmj,
                        self.v_coupled,
                        kind="cubic",
                        bounds_error=False,
                    )
                    igrid = np.linspace(
                        -self.npts_neg[self.inm] * self.fnmrms * self.nmrms[self.inm],
                        self.npts_pos[self.inm] * self.fnmrms * self.nmrms[self.inm],
                        self.nint,
                    )
                    jgrid = np.linspace(
                        -self.npts_neg[self.jnm] * self.fnmrms * self.nmrms[self.jnm],
                        self.npts_pos[self.jnm] * self.fnmrms * self.nmrms[self.jnm],
                        self.nint,
                    )
                    vijgrid = (
                        vtspl(igrid, jgrid)
                        - vtspl(igrid, jgrid * 0.0)
                        - vtspl(igrid * 0.0, jgrid)
                        + vtspl(igrid * 0.0, jgrid * 0.0)
                    )

                    # Save coupling correction to file for vistualisation.
                    info(
                        " @NM: Saving the interpolated potential energy to %s"
                        % (self.v_coupled_grid_filename,),
                        verbosity.medium,
                    )
                    outfile = self.imm.output_maker.get_output(
                        self.v_coupled_grid_filename, mode="wb"
                    )
                    np.save(outfile, vijgrid)
                    outfile.close_stream()

                    tmpfile = (
                        self.v_coupled_grid_file_prefix
                        + "."
                        + str(self.jnm)
                        + "."
                        + str(self.inm)
                        + ".dat"
                    )
                    # vtspl = interp2d(displacements_nmj, displacements_nmi, self.v_coupled, kind='cubic', bounds_error=False)
                    # vijgrid = vtspl(jgrid, igrid) - vtspl(jgrid, igrid * 0.0) - vtspl(jgrid * 0.0, igrid) + vtspl(jgrid * 0.0, igrid * 0.0)

                    # Save coupling correction to file for vistualisation.
                    info(
                        " @NM: Saving the interpolated potential energy to %s"
                        % (tmpfile,),
                        verbosity.medium,
                    )
                    outfile = self.imm.output_maker.get_output(tmpfile, mode="wb")
                    np.save(outfile, vijgrid.T)
                    outfile.close_stream()

        # Solves the SE once the mapping is finished.
        elif self.solve is True:
            self.solver()

        else:
            self.terminate()

    def solver(self):
        """
        Solves the VSCF equations in a mean-field manner.
        """

        # Saves the interpolated potential for each normal mode separately.
        for inm in self.inms:
            # Skips iteration if the file exists.
            ofn = self.v_coupled_grids_file_prefix + "." + str(inm) + ".dat"
            if os.path.exists(self.imm.output_maker.prefix + "." + ofn):
                continue

            # Stores the coupled potential
            vc = np.zeros((self.dof, self.nint, self.nint))
            for jnm in self.inms:
                if inm == jnm:
                    continue

                fn = (
                    self.v_coupled_grid_file_prefix
                    + "."
                    + str(inm)
                    + "."
                    + str(jnm)
                    + ".dat"
                )
                vc[jnm] = np.load(self.imm.output_maker.prefix + "." + fn).T

            info(
                " @NM: Saving the interpolated potentials for normal mode no. %8d in %s"
                % (inm, ofn),
                verbosity.medium,
            )
            outfile = self.imm.output_maker.get_output(ofn, mode="wb")
            np.save(outfile, vc)
            outfile.close_stream()

        # Initializes the independent mode free and internal energy.
        ai, ei = np.zeros(self.dof), np.zeros(self.dof)
        self.v_mft_grids = self.v_indep_grids.copy()

        info("\n", verbosity.medium)
        # Initializes the wavefunctions for all the normal modes.
        for inm in self.inms:
            for ibasis in range(self.nbasis):
                self.psi_i_grids[inm, ibasis, :] = self.psi(
                    ibasis, 1.0, self.imm.w[inm], self.q_grids[inm]
                )
                self.psi_i_grids[inm, ibasis, :] /= np.sqrt(
                    np.sum(self.psi_i_grids[inm, ibasis, :] ** 2)
                )

            (
                ai[inm],
                ei[inm],
                self.evals_imf[inm],
                self.evecs_imf[inm],
            ) = self.solve_schroedingers_equation(
                self.imm.w[inm], self.psi_i_grids[inm], self.v_indep_grids[inm], True
            )
            info(
                " @NM: The IMF free energy of mode %8d is %10.8e" % (inm, ai[inm]),
                verbosity.medium,
            )

        vscf_iter = 0
        a_imf = self.v0 + ai.sum()
        a_imf = ai.sum()
        a_vscf = a_imf
        self.evals_vscf, self.evecs_vscf = self.evals_imf.copy(), self.evecs_imf.copy()
        info(" @NM: The total IMF free energy is %10.8e" % (a_imf), verbosity.medium)
        info("\n @NM: The SCF begins.", verbosity.medium)

        while True:
            # Calculates the VSCF free energy.
            a_vscf_old = a_vscf
            a_vscf = self.v0 + ai.sum()
            a_vscf = ai.sum()
            vscf_iter += 1
            da = np.absolute(a_vscf - a_vscf_old) / len(self.inms)
            info(
                " @NM: CONVERGENCE : iteration = %8d   A =  %10.8e    D(A) = %10.8e / %10.8e"
                % (vscf_iter, a_vscf, da, self.athresh),
                verbosity.medium,
            )

            # Calculates the thermal density for each normal mode.
            # This is essentially the square of the wave function
            # times the Boltzmann density of each state.
            # info(' @NM: Calculating the thermal density.', verbosity.medium)
            self.rho_grids *= 0.0
            for inm in self.inms:
                for ibasis in range(self.nbasis):
                    self.rho_grids[inm] += (
                        np.exp(
                            -1.0
                            * (self.evals_vscf[inm, ibasis] - self.evals_vscf[inm, 0])
                            / self.imm.temp
                        )
                        * np.dot(self.psi_i_grids[inm].T, self.evecs_vscf[inm, ibasis])
                        ** 2
                    )
                    self.rho_grids[inm] /= self.rho_grids[inm].sum()

            # Calculates the mean field potential for each normal
            # mode and solves the Schroedinger's equation.
            for inm in self.inms:
                self.v_mft_grids[inm] = (
                    self.v_indep_grids[inm] + (1 - self.alpha) * self.v_mft_grids[inm]
                )
                fn = self.v_coupled_grids_file_prefix + "." + str(inm) + ".dat"
                vcg = np.load(self.imm.output_maker.prefix + "." + fn)
                for jnm in self.inms:
                    self.v_mft_grids[inm] += (
                        self.alpha * np.dot(vcg[jnm], self.rho_grids[jnm]) / self.nprim
                    )

                (
                    ai[inm],
                    ei[inm],
                    self.evals_vscf[inm],
                    self.evecs_vscf[inm],
                ) = self.solve_schroedingers_equation(
                    self.imm.w[inm], self.psi_i_grids[inm], self.v_mft_grids[inm], True
                )

            # Checks the convergence of the SCF procedure.
            da = np.absolute(a_vscf - a_vscf_old) / len(self.inms)
            if da < self.athresh and vscf_iter > 4:
                info("\n @NM: Convergence reached.", verbosity.medium)
                info(
                    " @NM: IMF free energy             = %10.8e" % (a_imf / self.nprim),
                    verbosity.low,
                )
                info(
                    " @NM: VSCF free energy correction = %10.8e"
                    % ((a_vscf - a_imf) / self.nprim),
                    verbosity.low,
                )
                info(
                    " @NM: ALL QUANTITIES PER PRIMITIVE UNIT CELL (WHERE APPLICABLE) \n",
                    verbosity.low,
                )
                info(
                    " @NM: Saving the energy eigenvalues in %s"
                    % (self.eigenvalue_filename),
                    verbosity.low,
                )
                outfile = self.imm.output_maker.get_output(self.eigenvalue_filename)
                np.savetxt(outfile, self.evals_vscf)
                outfile.close_stream()
                info(
                    " @NM: Saving the energy eigenvectors in %s in binary format.\n"
                    % (self.eigenvector_filename),
                    verbosity.low,
                )
                outfile = self.imm.output_maker.get_output(
                    self.eigenvector_filename, "wb"
                )
                np.save(outfile, self.evecs_vscf)
                outfile.close_stream()
                self.terminate()

    def one_dimensional_mapper(self, step):
        """
        Maps the potential energy landscape along a mode and returns
        the number of sampled points and the sampled potential energy.
        """

        # Determines sampling range for given normal mode
        nmd = self.fnmrms * self.imm.nmrms[self.inm]

        # Determines the displacement vector in Cartesian space.
        dev = np.real(self.imm.V.T[self.inm]) * nmd * np.sqrt(self.nprim)

        # Adds the minimum configuration to the list of sampled potential.
        v_indeps = [self.v0]

        # Displaces along the negative direction.
        counter = -1
        while True:
            self.imm.dbeads.q = self.imm.beads.q + dev * counter
            v = dstrip(self.imm.dforces.pots).copy()[0] / self.nprim
            v_indeps.append(v)

            # Bails out if the sampled potenital energy exceeds a user defined threshold.
            if self.nevib * self.imm.nmevib[self.inm] < np.absolute(v - self.v0):
                # Adds two extra points required later for solid spline fitting at edges.
                self.imm.dbeads.q -= dev
                v_indeps.append(dstrip(self.imm.dforces.pots).copy()[0] / self.nprim)
                self.imm.dbeads.q -= dev
                v_indeps.append(dstrip(self.imm.dforces.pots).copy()[0] / self.nprim)
                break

            counter -= 1

        r_npts_neg = -counter

        # invert vi so it starts with the potential for the most negative displacement and ends on the equilibrium
        v_indeps = v_indeps[::-1]

        counter = 1
        while True:
            self.imm.dbeads.q = self.imm.beads.q + dev * counter
            v = dstrip(self.imm.dforces.pots).copy()[0] / self.nprim
            v_indeps.append(v)

            # Bails out if the sampled potenital energy exceeds a user defined threshold.
            if self.nevib * self.imm.nmevib[self.inm] < np.absolute(v - self.v0):
                # add two extra points required later for solid spline fitting at edges
                self.imm.dbeads.q += dev
                v_indeps.append(dstrip(self.imm.dforces.pots).copy()[0] / self.nprim)
                self.imm.dbeads.q += dev
                v_indeps.append(dstrip(self.imm.dforces.pots).copy()[0] / self.nprim)
                break

            counter += 1

        r_npts_pos = counter

        return r_npts_neg, r_npts_pos, v_indeps

    def terminate(self):
        """
        Triggers a soft exit.
        """

        softexit.trigger(
            status="success", message=" @NM: The VSCF calculation has terminated."
        )
