#!/usr/bin/env python3

import numpy as np
import sys
import argparse
from ipi.utils.io import read_file
from ipi.utils.depend import *
from ipi.utils.units import *
from ipi.utils.messages import verbosity

description = """Computes the Transient Anisotropic Gaussian (TAG) approximation
of the instantaneous kinetic energy tensor, with a moving average
triangular window of the specified lag. Needs files with
the per-atom diagonal and off-diagonal components of the kinetic
energy tensor estimator.

Assumes the input files are in xyz format and atomic units,
with prefix.kin.xyz and prefix.kod.xyz naming scheme.
"""


def main(prefix, lag):
    verbosity.level = "low"
    lag = int(lag)

    # hard-coded prefixes for kinetic energy components  (kin=diagonal bit, kod=off-diagonal)
    ikin = open(prefix + ".kin.xyz", "r")
    ikod = open(prefix + ".kod.xyz", "r")
    otag = open(prefix + ".ktag_" + str(lag) + ".xyz", "w")

    natoms = 0
    ifr = 0
    cbuf = 2 * lag + 1
    while True:
        try:
            ret = read_file("xyz", ikin)
            tk = ret["atoms"]
            kin = dstrip(tk.q)
            ret = read_file("xyz", ikod)
            kod = dstrip(ret["atoms"].q)
            if natoms == 0:  # initializes vectors
                natoms = len(kin) / 3
                ktbuf = np.zeros(
                    (cbuf, natoms, 3, 3), float
                )  # implement the buffer as a circular one so one doesn't need to re-allocate and storage is continuous
                nkt = np.zeros((natoms, 3, 3), float)
                mkt = np.zeros((natoms, 3, 3), float)
                # akt = np.zeros((natoms, 3), float)
                mea = np.zeros((natoms, 3), float)
                mev = np.zeros((natoms, 3, 3), float)

            nkt[:, 0, 0] = kin[0 : natoms * 3 : 3]
            nkt[:, 1, 1] = kin[1 : natoms * 3 : 3]
            nkt[:, 2, 2] = kin[2 : natoms * 3 : 3]
            nkt[:, 0, 1] = nkt[:, 1, 0] = kod[0 : natoms * 3 : 3]
            nkt[:, 0, 2] = nkt[:, 2, 0] = kod[1 : natoms * 3 : 3]
            nkt[:, 2, 1] = nkt[:, 1, 2] = kod[2 : natoms * 3 : 3]
        except EOFError:  # Finished reading
            sys.exit(0)

        ktbuf[ifr % cbuf] = nkt
        if ifr >= (2 * lag):
            # now we can compute the mean tensor, and estimate the components of the kinetic energy
            mkt[:] = 0.0
            tw = 0.0
            for j in range(cbuf):
                w = 1.0 - np.abs(j - lag) * 1.0 / lag
                mkt += w * ktbuf[(ifr - j) % cbuf]
                tw += w
            mkt *= 1.0 / tw

            # computes the eigenvalues of the TAG tensor
            for i in range(natoms):
                [mea[i], mev[i]] = np.linalg.eigh(mkt[i])

            otag.write(
                "%d\n# TAG eigenvalues e1 e2 e3 with lag %d. Frame: %d\n"
                % (natoms, lag, ifr - lag)
            )
            for i in range(natoms):
                otag.write(
                    "%6s  %15.7e  %15.7e  %15.7e\n"
                    % (tk.names[i], mea[i, 0], mea[i, 1], mea[i, 2])
                )

        ifr += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "prefix", type=str, nargs=1, help="Prefix of the trajectories to process."
    )
    parser.add_argument(
        "lag", type=int, nargs=1, help="Time lag for the TAG, in number of time steps."
    )

    args = parser.parse_args()

    main(args.prefix[0], args.lag[0])
