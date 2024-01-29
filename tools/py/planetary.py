#!/usr/bin/env python3

import sys
import argparse
import importlib

import numpy as np
from numpy import matmul
from numpy.linalg import multi_dot

from ipi.utils import sparse
from ipi.utils.io import read_file
from ipi.utils.prng import Random
from ipi.utils.messages import warning
from ipi.utils.depend import *
from ipi.utils.units import *
from ipi.utils.io import netstring_encoded_loadz

# *******************************#

description = """\
---------------
Planetary model
---------------

Process data from a single centroid trajectory to obtain TCFs using the
planetary model.

Assumes existance of the following files in the working directory:

    <prefix>.xc.xyz
  Stores centroid configurations along trajectory

    <prefix>.pc.xyz
  Stores centroid momenta along trajectory

    <prefix>.omega2
  Binary file storing path-integral frequency matrices at each centroid
  configuration along trajectory

User must provide a list of external modules as a command line argument,
which the program will use to work out how to calculate the TCFs. A
template for constructing a module in the required format is given in
ipi/tools/py/estmod_example.py.

Written by Raz Benson

"""


def getInput():
    """
    Parse command line arguments and return input parameters for Planets class.

    Returns
        prefix : filename prefix
        natoms : number of atoms
        iseed : random number seed
        npl : number of planets to simulate
        npts : number of data points
        stride : number of planetary steps between centroid updates
        dt : planetary time step in fs
        temperature : temperature in K
        estimators : list of modules defining estimators in the required format

    """

    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("prefix", type=str, help="Filename prefix")

    parser.add_argument(
        "dt",
        type=float,
        help="Time step (fs) for planets. Does NOT have to be the same as for the centroids.",
    )

    parser.add_argument(
        "temperature",
        type=float,
        help="Simulation temperature (K). Best to make sure it is the same as the centroid temperature!",
    )

    parser.add_argument(
        "stride",
        type=int,
        help="""Number of steps to take between centroid updates
                             (so stride*dt should be the time elapsed between centroid positions being RECORDED).""",
    )

    parser.add_argument("npl", type=int, help="Number of planets for each centroid")

    parser.add_argument(
        "npts", type=int, help="Number of points in centroid trajectory"
    )

    parser.add_argument("iseed", type=int, help="Random number seed")

    parser.add_argument(
        "estmods",
        type=str,
        nargs="+",
        help="External modules defining estimators in the required format",
    )

    args = parser.parse_args()

    with open("{}.xc.xyz".format(args.prefix), "r") as f:
        args.natoms = int(f.readline())
        for i in range(args.natoms + 2):
            f.readline()
        stride2 = int(f.readline().split("Step:")[1].split("Bead:")[0])

    if args.stride != stride2:
        warning(
            "Make sure stride*dt = time elapsed between centroid positions being recorded\n"
            + "(this might already be true)."
        )

    estimators = []
    for modname in args.estmods:
        # This might be quite ugly but it works for now
        if "/" in modname:
            path, modname = modname.rsplit("/", 1)
            sys.path.insert(0, path)
        mod = importlib.import_module(modname)
        # estimators.append(mod.est)
        estimators.append(mod)

    args.estmods = ", ".join(args.estmods)
    print("----------------------")
    print("DETERMINED INPUT PARAMETERS")
    print("----------------------")
    for k, v in sorted(args.__dict__.items()):
        print("{:<30s}{}".format(k, v))
    print("----------------------")
    sys.stdout.flush()

    return (
        args.prefix,
        args.natoms,
        args.iseed,
        args.npl,
        args.npts,
        args.stride,
        args.dt,
        args.temperature,
        estimators,
    )


# *******************************#


def simple_corr(A, B, **kwargs):
    """
    Simplest correlater applicable to vector opeators e.g. dipole autocorrelation function

    Args:
        A : array with shape (npts, npl, ndim)
        B : array with shape (npts, npl, ndim)
    Returns:
        tcf : 1-d array of length npts

    Here npts is number of data points in trajectory, npl is number planets and ndim is
    number of dimensions (e.g. ndim = 3 if A and B represent dipole moment operators).
    """
    npts, npl, ndim = A.shape
    tcf = np.zeros(npts)
    weights = np.arange(npts, 0, -1, dtype=float)
    for i in range(npl):
        for j in range(ndim):
            tcf[:] += (
                np.correlate(A[:, i, j], B[:, i, j], mode="full")[npts - 1 :: -1]
                / weights
            )
    # for i in xrange(npl):
    #    for j in xrange(ndim):
    #        tcf[0] += np.mean(A[:,i,j] * B[:,i,j])
    #        for k in xrange(1, npts):
    #            tcf[k] += np.mean(A[0:-k,i,j] * B[k:,i,j])
    tcf /= float(npl)
    return tcf


# *******************************#


class Planets(object):
    """
    Stores relevant trajectory info for planetary model simulation, calculates
    planetary dynamics and evaluates estimators and TCFs.
    """

    def __init__(
        self, prefix, natoms, iseed, npl, npts, stride, dt, temperature, estimators
    ):
        """
        Initialises planetary model simulation.

        Args:
            prefix : identifies trajectory files
            natoms : number of atoms
            iseed : random number seed
            stride : number of planetary steps between centroid updates
            dt : planetary time step in fs
            temperature : temperature in K
            estimators : list of modules defining estimators in the required format
                         (see ipi/tools/py/estmod_example.py)

        """
        self.prefix = prefix
        self.natoms = natoms
        self.iseed = iseed
        self.npl = npl
        self.npts = npts
        self.stride = stride
        self.dt = unit_to_internal("time", "femtosecond", dt)
        self.temperature = unit_to_internal("temperature", "kelvin", temperature)
        self.estimators = estimators
        self.prng = Random(self.iseed)
        self.beta = 1.0 / self.temperature

        # Centroid
        self.qc = np.zeros(3 * self.natoms)
        self.pc = np.zeros(3 * self.natoms)
        # Fluctuation (planet - centroid)
        self.p = np.zeros((3 * self.natoms, self.npl))
        self.q = np.zeros((3 * self.natoms, self.npl))
        # Dimensionless fluctuation
        self.qtil = np.zeros((3 * self.natoms, self.npl))
        self.ptil = np.zeros((3 * self.natoms, self.npl))
        # Total (planet coorindates)
        self.psum = np.zeros((3 * self.natoms, self.npl))
        self.qsum = np.zeros((3 * self.natoms, self.npl))
        # Path-integral frequency
        self.omega2 = np.zeros((3 * self.natoms, 3 * self.natoms))
        self.omega = np.zeros((3 * self.natoms, 3 * self.natoms))
        self.omega_old = np.zeros((3 * self.natoms, 3 * self.natoms))
        self.omega_interp = np.zeros((3 * self.natoms, 3 * self.natoms))
        # omega2 eigensystem
        self.evals = np.zeros(3 * self.natoms)
        self.evals_sqrt = np.zeros(3 * self.natoms)
        self.evecs = np.zeros((3 * self.natoms, 3 * self.natoms))
        # Position smearing matrix
        self.a = np.zeros((3 * self.natoms, 3 * self.natoms))
        self.a_sqrt = np.zeros((3 * self.natoms, 3 * self.natoms))
        self.a_inv = np.zeros((3 * self.natoms, 3 * self.natoms))
        # Momentum smearing matrix
        self.b = np.zeros((3 * self.natoms, 3 * self.natoms))
        self.b_sqrt = np.zeros((3 * self.natoms, 3 * self.natoms))
        self.b_inv = np.zeros((3 * self.natoms, 3 * self.natoms))
        # Boolean mask for frequencies close to 0
        self.mask = np.zeros(3 * self.natoms, dtype=bool)

        self.fomega2 = open("{}.omega2".format(self.prefix), "r")
        self.fqc = open("{}.xc.xyz".format(self.prefix), "r")
        self.fpc = open("{}.pc.xyz".format(self.prefix), "r")

        for i, mod in enumerate(self.estimators):
            if not hasattr(mod, "corr"):
                mod.corr = simple_corr
            mod.TCF = np.zeros(self.npts)
            mod.TCF_c = np.zeros(self.npts)
            mod.Aarr = np.zeros((self.npts, self.npl) + mod.Ashape)
            mod.Aarr_c = np.zeros((self.npts, 1) + mod.Ashape)
            mod.Atemp = np.zeros(np.prod(mod.Ashape))
            mod.Barr = np.zeros((self.npts, self.npl) + mod.Bshape)
            mod.Barr_c = np.zeros((self.npts, 1) + mod.Bshape)

    def sample(self):
        """
        Sample planetary momenta in dimensionless normal mode coordinates
        """
        self.qtil[:] = self.prng.gvec((3 * self.natoms, self.npl))
        self.ptil[:] = self.prng.gvec((3 * self.natoms, self.npl))

    def read_omega2(self):
        """
        Read in next instance of path-integral frequency matrix from <prefix>.omega2
        """
        omega2 = sparse.load_npz(self.fomega2, loader=netstring_encoded_loadz)
        # Convert self.omega2 to double precision (need at least single precision for np.linalg to work)
        self.omega2[:] = omega2.toarray().astype(np.float64)

    def read_qcpc(self):
        """
        Read in next instances of centroid positions and momenta
        """
        ret = read_file("xyz", self.fqc)  # , readcell=True)
        self.qpos = ret["atoms"]
        ret = read_file("xyz", self.fpc)  # , readcell=True)
        self.ppos = ret["atoms"]
        self.qc[:] = self.qpos.q
        self.pc[:] = self.ppos.q

    def matrix_setup(self):
        """
        Set up and store as attributes the various matrices that need to be manipulated
        """
        self.omega[:] = 0.0
        self.a[:] = 0.0
        self.a_inv[:] = 0.0
        self.a_sqrt[:] = 0.0
        self.b[:] = 0.0
        self.b_inv[:] = 0.0
        self.b_sqrt[:] = 0.0

        self.evals[:], self.evecs[:] = np.linalg.eigh(self.omega2, UPLO="L")

        self.mask[:] = self.beta**2 * self.evals < 1e-14
        self.evals[self.mask] = 0.0
        self.evals_sqrt[:] = np.sqrt(self.evals)
        np.fill_diagonal(self.omega, self.evals_sqrt)

        # Calculate approximate smearing matrices
        self.a[~self.mask, ~self.mask] = (
            0.5
            * self.beta
            * self.evals_sqrt[~self.mask]
            / np.tanh(0.5 * self.beta * self.evals_sqrt[~self.mask])  # <=> *coth()
            - 1.0
        ) / (self.beta * self.evals[~self.mask])
        self.a[self.mask, self.mask] = self.beta / 12.0
        np.fill_diagonal(self.a_inv, 1.0 / np.diag(self.a))
        np.fill_diagonal(self.a_sqrt, np.sqrt(np.diag(self.a)))

        self.b[~self.mask, ~self.mask] = (
            self.evals[~self.mask] * self.a[~self.mask, ~self.mask]
        )  # *self.m[~self.mask]**2
        self.b_inv[~self.mask, ~self.mask] = 1.0 / self.b[~self.mask, ~self.mask]
        np.fill_diagonal(self.b_sqrt, np.sqrt(np.diag(self.b)))

        # Undiagonalize matrices inplace
        self.eigenbasis(self.omega, backward=True)
        self.eigenbasis(self.a, backward=True)
        self.eigenbasis(self.a_sqrt, backward=True)
        self.eigenbasis(self.a_inv, backward=True)
        self.eigenbasis(self.b, backward=True)
        self.eigenbasis(self.b_sqrt, backward=True)
        self.eigenbasis(self.b_inv, backward=True)

        # Ensure appropriate scaling by mass matrix - very important
        # to to this in the correct basis!
        self.a[:] = multi_dot([self.mmat_sqrt_inv, self.a, self.mmat_sqrt_inv])
        self.a_sqrt[:] = multi_dot([self.mmat_qrt_inv, self.a_sqrt, self.mmat_qrt_inv])
        self.a_inv[:] = multi_dot([self.mmat_sqrt, self.a_inv, self.mmat_sqrt])
        self.b[:] = multi_dot([self.mmat_sqrt, self.b, self.mmat_sqrt])
        self.b_sqrt[:] = multi_dot([self.mmat_qrt, self.b_sqrt, self.mmat_qrt])
        self.b_inv[:] = multi_dot([self.mmat_sqrt_inv, self.b_inv, self.mmat_sqrt_inv])

    def eigenbasis(self, mat, backward=False):
        """
        Transfrom a given matrix in place to or from the eigenbasis of the
        path-integral frequency matrix.

        Args:
            mat : matrix with shape (3*natoms, 3*natoms)
            backward : boolean specifying whether to transform to or from eigenbasis

        """
        if backward:
            mat[:] = multi_dot([self.evecs, mat, self.evecs.transpose()])
        else:
            mat[:] = multi_dot([self.evecs.transpose(), mat, self.evecs])

    def step(self):
        """
        Evolve planets for time stride*dt
        """
        self.omega_old[:] = self.omega.copy()
        self.read_omega2()
        self.read_qcpc()
        self.matrix_setup()

        for j in range(self.stride):
            # Linear interpolation of frequency matrix
            self.omega_interp[:] = (
                j * self.omega + (self.stride - j) * self.omega_old
            ) / self.stride
            # Velocity verlet integrator
            self.ptil[:] -= 0.5 * self.dt * matmul(self.omega_interp, self.qtil)
            self.qtil[:] += self.dt * matmul(self.omega_interp, self.ptil)
            self.ptil[:] -= 0.5 * self.dt * matmul(self.omega_interp, self.qtil)

        self.q[:] = multi_dot(
            [self.mmat_qrt_inv, self.a_sqrt, self.mmat_qrt, self.qtil]
        )
        self.p[:] = multi_dot(
            [self.mmat_qrt, self.b_sqrt, self.mmat_qrt_inv, self.ptil]
        )
        self.qsum[:] = self.q + self.qc[:, np.newaxis]
        self.psum[:] = self.p + self.pc[:, np.newaxis]

    def estimate(self, i):
        """
        Evaluate TCF estimators.

        Args:
            i : index specifying number of centroid samples elapsed

        """
        for mod in self.estimators:
            if hasattr(mod, "Afunc0"):
                # Might as well calculate centroid TCF at little extra cost
                mod.Barr_c[i, :] += mod.Bfunc(
                    self.qc[:, np.newaxis], self.pc[:, np.newaxis]
                )
                mod.Aarr_c[i, :] += mod.Afunc0(
                    self.qc[:, np.newaxis], self.pc[:, np.newaxis]
                )

            # Deal with B(Q_0+q) first, this is easy
            mod.Barr[i, :] += mod.Bfunc(self.qsum, self.psum)

            # Now deal with f_A(Q0+q,p)

            if mod.method == "0thorder_re":
                # 0th order term (evaluate the function A at the planet position)
                mod.Aarr[i, :] += mod.Afunc0(self.qsum, self.psum)

            elif mod.method == "1storder_im":
                # 1st order term, with shape given by
                # ( npl,    prod(mod.Ashape),    3*self.natoms)
                A1q = mod.Afunc1(self.qsum)

                for j in range(self.npl):
                    mod.Atemp[:] = 0.0

                    for k, val in enumerate(A1q[j]):
                        mod.Atemp[k] = multi_dot([self.p[:, j], self.b_inv, val])

                    mod.Aarr[i, j, :] += mod.Atemp.reshape(mod.Ashape) / 2.0

            elif mod.method == "2ndorder_re":
                # 0th order term (evaluate the function A at the planet position)
                mod.Aarr[i, :] += mod.Afunc0(self.qsum, self.psum)
                # 2nd order term, with shape given by
                # ( npl,    prod(mod.Ashape),    3*self.natoms,    3*self.natoms )
                A2q = mod.Afunc2(self.qsum)

                for j in range(self.npl):
                    mod.Atemp[:] = 0.0

                    for k, val in enumerate(A2q[j]):
                        mod.Atemp[k] = np.trace(matmul(self.b_inv, val)) - multi_dot(
                            [self.p[:, j], self.b_inv, val, self.b_inv, self.p[:, j]]
                        )

                    mod.Aarr[i, j, :] += mod.Atemp.reshape(mod.Ashape) / 8.0

            else:
                raise ValueError(
                    "f_A approximation method '{}' not recognised".format(mod.method)
                )

    def correlate(self):
        """
        Run at the end to calculate the TCFs
        """
        for mod in self.estimators:
            mod.TCF[:] = mod.corr(mod.Aarr, mod.Barr, **self.__dict__)
            mod.TCF_c[:] = mod.corr(mod.Aarr_c, mod.Barr_c, **self.__dict__)

    def write_tcfs(self):
        """
        Run at the end to write the TCFs to .dat files
        """
        tarr = np.arange(self.npts) * self.dt * self.stride
        for mod in self.estimators:
            if hasattr(mod, "Afunc0"):
                np.savetxt(
                    "{}_ce_{}.dat".format(self.prefix, mod.name),
                    np.transpose([tarr, mod.TCF_c]),
                )
                print(
                    "Saved {} (centroid) to {}_ce_{}.dat".format(
                        mod.name, self.prefix, mod.name
                    )
                )
            np.savetxt(
                "{}_pl_{}.dat".format(self.prefix, mod.name),
                np.transpose([tarr, mod.TCF]),
            )
            print(
                "Saved {} (full) to {}_pl_{}.dat".format(
                    mod.name, self.prefix, mod.name
                )
            )

    def shutdown(self):
        """
        Prepare to end simulation by safely closing any open files and removing temporary files
        """
        # os.remove("TEMP2_PLANETARY")
        self.fomega2.close()
        self.fqc.close()
        self.fpc.close()

    def get_masses(self):
        """
        Deterine mass matrix based on identities of atoms
        """
        m = dstrip(self.qpos.m)
        self.m = np.concatenate((m, m, m))
        self.m[:] = self.m.reshape(3, -1).transpose().reshape(-1)
        self.mmat = np.diag(self.m)
        self.mmat_sqrt = np.sqrt(self.mmat)
        self.mmat_qrt = np.sqrt(self.mmat_sqrt)
        self.mmat_inv = np.linalg.inv(self.mmat)
        self.mmat_sqrt_inv = np.linalg.inv(self.mmat_sqrt)
        self.mmat_qrt_inv = np.linalg.inv(self.mmat_qrt)

    def simulation(self, correlate=True, write=True, shutdown=True):
        """
        Execute this method to run the planetary simulation.

        Args:
            correlate : set to True to evaluate TCFs after full time elapsed
            write : set to True to write TCFs to .dat files
            shutdown : set to True to safely shutdown after simulation complete
                       by calling self.shutdown()

        """
        print("STARTING PLANETARY SIMULATION")
        self.read_omega2()
        self.read_qcpc()
        self.get_masses()
        self.matrix_setup()
        self.sample()
        self.q[:] = multi_dot(
            [self.mmat_qrt_inv, self.a_sqrt, self.mmat_qrt, self.qtil]
        )
        self.p[:] = multi_dot(
            [self.mmat_qrt, self.b_sqrt, self.mmat_qrt_inv, self.ptil]
        )
        self.qsum[:] = self.q + self.qc[:, np.newaxis]
        self.psum[:] = self.p + self.pc[:, np.newaxis]
        self.estimate(0)
        self.count = 0
        for i in range(self.npts - 1):
            self.step()
            self.count += 1  # this might be pointless, obviously self.count = i
            self.estimate(i + 1)
            print("Completed step {:d} of {:d}".format(self.count, self.npts - 1))
            sys.stdout.flush()
        if correlate:
            print("CALCULATING TCFS")
            self.correlate()
        if write:
            print("SAVING TCFS")
            self.write_tcfs()
        if shutdown:
            print("SIMULATION COMPLETE. SHUTTING DOWN.")
            self.shutdown()
        else:
            warning("Reached end of simulation but files may remain open.")


def main():
    P = Planets(*getInput())
    P.simulation()


if __name__ == "__main__":
    main()
