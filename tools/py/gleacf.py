#!/usr/bin/env python3
"""
A command-line utility to apply the GLE/Langevin deconvolution scheme from
Rossi, Kapil, Ceriotti, JCP (2018) http://doi.org/10.1063/1.4990536
"""

import sys
import argparse
import numpy as np
from ipi.utils.tools import gle_frequency_kernel, isra_deconvolute, get_gle_matrices
from ipi.utils.messages import info, verbosity


def input_facf(path2inputfile, mrows, stride):
    """Reads the facf file."""
    dfacf = np.genfromtxt(path2inputfile, usecols=((0, 1)))
    if mrows == -1:
        mrows = len(dfacf)
    return dfacf[:mrows][::stride]


def output_facf(xy, oprefix, reffacf):
    """Exports the facf file."""
    np.savetxt(oprefix + "_convacf.data", np.vstack((xy[0], xy[1])).T)


def gleacf(
    path2ixml,
    path2iA,
    path2iC,
    path2ifacf,
    oprefix,
    action,
    nrows,
    stride,
    tscale,
    dparam,
):
    if path2ixml is not None and path2iA is not None:
        raise Exception(
            "The drift and covariance matrices have to be provided either through the i-pi xml file or manually. Can not use both forms of input simultaneously. "
        )
    elif path2ixml is not None and path2iA is None:
        Ap, Cp, Dp = get_gle_matrices(path2ixml, tscale)
    elif path2ixml is None and path2iA is not None:
        Ap = np.loadtxt(path2iA, dtype=float, ndmin=2) * tscale
        if path2iC is not None:
            Cp = np.loadtxt(path2iC)
        else:
            Cp = np.eye(len(Ap))
        Dp = np.dot(Ap, Cp) + np.dot(Cp, Ap.T)
    else:
        raise Exception(
            "The drift and covariance matrixes have to be provided either through the i-pi xml file or manually."
        )

    # imports the facf function.
    ifacf = input_facf(path2ifacf, nrows, stride)
    ix = ifacf[:, 0]
    iy = ifacf[:, 1]

    # computes the facf kernel
    info("# computing the kernel.", verbosity.low)
    ker = gle_frequency_kernel(ix, Ap, Dp)
    print(ker)

    # (de-)convolutes the spectrum
    if action == "conv":
        info("# printing the output spectrum.", verbosity.low)
        # yes, the convolution is just np.dot(iy, ker.T)
        output_facf((ix, np.dot(iy, ker.T)), oprefix, input_facf(path2ifacf, nrows, 1))
    elif action == "deconv":
        info("# deconvoluting the input spectrum.", verbosity.low)
        dey, ldey, err, lap = isra_deconvolute(ix, iy, ker, dparam[0], dparam[1])
        np.savetxt(oprefix + "_decacf.data", np.vstack([ix, dey]).T)
        np.savetxt(oprefix + ".err_lap", np.vstack([err, lap]).T)
        print(ldey.shape)
        np.savetxt(oprefix + ".history", ldey.T)


if __name__ == "__main__":
    # adds description of the program.
    parser = argparse.ArgumentParser(
        description="""Given the parameters of a Generalized Langevin Equation and the microcanonical 
                    density of states predicts the power spectrum of the velocity-velocity auto-correlation obtained by the dynamics. 
                    Conversely, given the power spectrum velocity-velocity autocorrelation function removes the disturbance affected 
                    by the thermostat and returns the underlying vibrational density of states. """
    )

    # adds arguments.
    parser.add_argument(
        "-a",
        "--action",
        choices=["conv", "deconv"],
        default=None,
        required=True,
        help="choose conv if you want to obtain the response of the thermostat on the vibrational density of states; choose deconv if you want to obtain the micro-canonical density of states by removing the disturbance induced by the thermostat",
    )
    parser.add_argument(
        "-ifacf",
        "--input_facf",
        type=str,
        default=None,
        help="The input Fourier transform of the auto-correlation function (FACF). The unit of time is assumed to be atomic units unless the additional -time_scaling parameter is provided.",
    )
    parser.add_argument(
        "-ixml",
        "--input_xml",
        type=str,
        default=None,
        help="The i-PI xml file that contains the Drift (A) and the Covariance (C) matrices as well as the ensemble temperature (T). If a normal mode GLE thermostat is used then only the parameters of the thermostat that acts on the centroid is are used assuming that the spectrum that is to be deconvoluted is computed from the dynamics of the centroid.  The units are interpretted from the file and the A and C matrices are converted to atomic units. The C matrix is further divided by k_B T. If the -time_scaling parameter is provided, the units of A are modified using the scaling factor.",
    )
    parser.add_argument(
        "-ia",
        "--input_A",
        type=str,
        default=None,
        help="The file containing the Drift matrix A of the GLE thermostat. The units are assumed to be atomic units. If the -time_scaling factor is provided the units of A are changed to match those of the FACF.",
    )
    parser.add_argument(
        "-ic",
        "--input_C",
        type=str,
        default=None,
        help="The file containing the dimensionless Covariance matrix C which is assumed to be rescaled with respect to k_B T. It is assumed to be an Identity matrix of the appropriate rank by default.",
    )
    parser.add_argument(
        "-mr",
        "--max_rows",
        type=int,
        default=-1,
        help="The index of the last row to be imported from FACF. Allows one to choose the maximum frequency which is considered during the process of (de)convolution.",
    )
    parser.add_argument(
        "-s",
        "--stride",
        type=int,
        default=1,
        help="The stride for importing the FACF. Allows one to choose the granularity of the frequency scale. ",
    )
    parser.add_argument(
        "-ts",
        "--time_scaling",
        type=float,
        default="1.0",
        help="The unit of time associated with the FACF w.r.t atomic units. Defaults to 1.0. ",
    )
    parser.add_argument(
        "-dp",
        "--deconv_param",
        nargs=2,
        type=int,
        default=[500, 10],
        help="The parameters associated with the deconvolution. Since the operation is based on an iterative algorithm, it requires the total number of epochs nepocs and the stride pstride at which the output spectrum is returned. Usage: -dp nepochs pstride",
    )
    parser.add_argument(
        "-op",
        "--output_prefix",
        type=str,
        default="output",
        help="the prefix of the (various) output files.",
    )

    # parses arguments.
    if len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        parser.print_help()
        sys.exit()

    # stores the arguments
    path2ixml = args.input_xml
    path2iA = args.input_A
    path2iC = args.input_C
    path2ifacf = args.input_facf
    oprefix = args.output_prefix
    action = args.action
    nrows = args.max_rows
    stride = args.stride
    tscale = args.time_scaling
    dparam = np.asarray(args.deconv_param, dtype=int)

    gleacf(
        path2ixml,
        path2iA,
        path2iC,
        path2ifacf,
        oprefix,
        action,
        nrows,
        stride,
        tscale,
        dparam,
    )
