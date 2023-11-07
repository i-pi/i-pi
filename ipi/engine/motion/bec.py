__all__ = ["BECTensorsCalculator"]

import numpy as np

from ipi.engine.motion import Motion
from ipi.utils.depend import *
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, info, warning
from ipi.utils.units import Constants


class BECTensorsCalculator(Motion):

    """BEC tensors calculator using finite difference."""

    def __init__(
        self,
        atoms=["all"],
        pos_shift=0.001,
        bec=np.zeros(0, float),
        prefix="",
        asr="none",
    ):
        """Initialises BECTensorsCalculator."""

        super(BECTensorsCalculator, self).__init__(fixcom=False, fixatoms=None)

        # self.mode = mode
        self.phcalc = FDBECTensorsCalculator()

        self.deltax = pos_shift
        self.bec = bec.copy()
        self.correction = np.zeros(0, float)

        self.prefix = prefix
        self.asr = asr
        self.atoms = atoms

        if self.prefix == "":
            self.prefix = "BEC"

    def bind(self, ens, beads, nm, cell, bforce, prng, omaker):
        super(BECTensorsCalculator, self).bind(
            ens, beads, nm, cell, bforce, prng, omaker
        )

        # Raises error for nbeads not equal to 1.
        if self.beads.nbeads > 1:
            raise ValueError(
                "Calculation not possible for number of beads greater than one."
            )

        self.tomove = np.full(3 * self.beads[0].natoms, False)
        if len(self.atoms) == 0:
            self.tomove = np.full(3 * self.beads[0].natoms, False)
        elif len(self.atoms) == 1 and self.atoms[0].lower() == "all":
            self.tomove = np.full(3 * self.beads[0].natoms, True)
        else:
            for i in self.atoms:
                if i.isdigit():
                    i = int(i)
                    self.tomove[i * 3 : (i + 1) * 3] = True
                else:
                    if i not in list(self.beads[0].names):
                        raise ValueError("wrong input")
                    else:
                        index = list(self.beads[0].names).index(i)
                    if not hasattr(index, "__len__"):
                        index = [index]
                    for j in index:
                        j = int(j)
                        self.tomove[j * 3 : (j + 1) * 3] = True

        self.phcalc.bind(self)

    def step(self, step=None):
        """Executes one step of BEC computation."""
        if step < 3 * self.beads.natoms:
            self.phcalc.step(step)
        else:
            self.phcalc.transform()
            self.apply_asr()
            self.printall()
            softexit.trigger(
                status="success",
                message="BEC tensors have been calculated. Exiting simulation",
            )

    def printall(self):
        """Prints matrices to file"""

        file = "{:s}.txt".format(
            self.prefix,
        )
        np.savetxt(file, self.bec.reshape((-1, 3)), delimiter=" ", fmt="%15.10f")

        if self.correction.shape != (0,):
            file = "{:s}.correction.txt".format(self.prefix)
            np.savetxt(
                file, self.correction.reshape((-1, 3)), delimiter=" ", fmt="%15.10f"
            )

        return

    def apply_asr(self):
        """
        Removes the translations and/or rotations depending on the asr mode.
        """

        #
        # Translational Sum Rule :
        # \sum_I Z^I_ij = 0
        #
        # This means that the translation of all the ions does not lead to any change in the dipole.
        #
        # So we compute this sum, which should be zero, but it is not due to "numerical noise",
        # and then we subtract this amount (divided by the number of atoms) to each BEC.
        #
        # Pay attention that in this case self.bec has already the shape (natoms,3,3)
        # and self.correction has the shape (3,3)
        #

        if self.asr == "lin":
            if np.all(self.tomove):
                warning("Sum Ruls can not be applied because some dofs were kept fixed")

            self.correction = self.bec.sum(axis=0) / self.beads.natoms
            self.bec -= self.correction

        elif self.asr == "none":
            return

        # We should add the Rotational Sum Rule(s)


class FDBECTensorsCalculator(dobject):

    """Finite difference BEC tensors evaluator."""

    def __init__(self):
        pass

    def bind(self, dm):
        """Reference all the variables for simpler access."""

        self.dm = dm
        self.original = np.asarray(dstrip(self.dm.beads.q[0]).copy())
        self.atoms = dm.atoms

        def check_dimension(M):
            if M.size != (3 * self.dm.beads.q.size):
                if M.size == 0:
                    M = np.full((self.dm.beads.q.size, 3), np.nan, dtype=float)
                else:
                    raise ValueError("matrix size does not match system size")
            else:
                M = M.reshape(((self.dm.beads.q.size, 3)))
            return M

        # Initialises a 3*number of atoms X 3*number of atoms dynamic matrix.
        self.dm.bec = check_dimension(self.dm.bec)

        return

    def step(self, step=None):
        """Computes one row of the BEC tensors"""

        if self.dm.tomove[step]:
            # initializes the finite deviation
            dev = np.zeros(3 * self.dm.beads.natoms, float)

            # displacement in cartesian components
            dev[step] = self.dm.deltax

            # displaces kth d.o.f by delta.
            self.dm.beads.q.set(self.original + dev)
            # Tplus = np.asarray(dstrip(self.dm.ensemble.eda.polarization).copy())

            # get the dipole, which is defined also for isolated systems
            Dplus = np.asarray(dstrip(self.dm.ensemble.eda.dipole).copy())

            # displaces kth d.o.f by -delta.
            self.dm.beads.q.set(self.original - dev)
            # Tminus = np.asarray(dstrip(self.dm.ensemble.eda.polarization).copy())

            # get the dipole, which is defined also for isolated systems
            Dminus = np.asarray(dstrip(self.dm.ensemble.eda.dipole).copy())

            # self.dm.bec[step] = ( self.dm.ensemble.cell.V / Constants.e ) * ( Tplus - Tminus ) / ( 2 * self.dm.deltax )
            self.dm.bec[step] = (
                (1.0 / Constants.e) * (Dplus - Dminus) / (2 * self.dm.deltax)
            )

        else:
            info(" We have skipped the dof # {}.".format(step), verbosity.low)

        return

    def transform(self):
        # reshape
        self.dm.bec = self.dm.bec.reshape((self.dm.beads.natoms, 3, 3))

        return
