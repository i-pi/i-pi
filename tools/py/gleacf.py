#!/usr/bin/env python3
"""

Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.

Cuts short the output of a previous i-pi simulation, up to the
step indicated in the <step> field of the input file.
This is useful to restart a simulation that crashed.

It should be run in the same dyrectory as where i-pi was (or is being)
run, and simply fetches all information from the simulation input file.
One should also specify a directory name in which the trimmed files
will be output.

Syntax:
   trimsim.py inputfile.xml
"""


import sys
import argparse
import numpy as np
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml


def input_facf(path2inputfile, mrows, stride):
    """Reads the facf file."""
    dfacf = np.genfromtxt(path2inputfile, usecols=((0, 1)))
    if mrows == -1:
        mrows = len(dfacf)
    return dfacf[:mrows][::stride]


def output_facf(xy, oprefix, reffacf):
    """Exports the facf file."""
    np.savetxt(oprefix + "_facf.data", np.vstack((xy[0], xy[1])).T)


def Aqp(omega_0, Ap):
    """Given the free particle Ap matrix and the frequency of the harmonic oscillator, computes the full drift matrix."""
    dAqp = np.zeros(np.asarray(Ap.shape) + 1)
    dAqp[1, 0] = -np.power(omega_0, 2)
    dAqp[0, 1] = 1
    dAqp[1:, 1:] = Ap
    return dAqp


def Dqp(omega_0, Dp):
    """Given the free particle Dp matrix and the frequency of the harmonic oscillator, computes the full D matrix."""
    dDqp = np.zeros(np.asarray(Dp.shape) + 1)
    dDqp[1:, 1:] = Dp
    return dDqp


# def Cqp(omega_0, Ap, Dp):
#    """Given the free particle Ap and Dp matrices and the frequency of the harmonic oscillator, computes the full covariance matrix."""
#    dAqp = Aqp(omega_0, Ap)
#    dDqp = Dqp(omega_0, Dp)
#    return sp.solve_continuous_are(
#        -dAqp, np.zeros(dAqp.shape), dDqp, np.eye(dAqp.shape[-1])
#    )


def Cqp(omega0, idAqp, idDqp):
    """Given the free particle Ap and Dp matrices and the frequency of the harmonic oscillator, computes the full covariance matrix."""
    # "stabilizes" the calculation by removing the trivial dependence of <a^2> on omega0 until the very end
    dAqp = idAqp.copy()
    dDqp = idDqp.copy()
    dAqp[:, 0] *= omega0
    dAqp[0, :] /= omega0
    dDqp[:, 0] /= omega0
    dDqp[:, 0] /= omega0

    # solve "a la MC thesis" using just numpy
    a, omat = np.linalg.eig(dAqp)
    oinv = np.linalg.inv(omat)
    W = np.dot(np.dot(oinv, dDqp), oinv.T)
    for i in range(len(W)):
        for j in range(len(W)):
            W[i, j] /= a[i] + a[j]
    nC = np.real(np.dot(omat, np.dot(W, omat.T)))

    nC[:, 0] /= omega0
    nC[0, :] /= omega0
    return nC


def gleKernel(omega, Ap, Dp):
    """Given the Cp and Dp matrices for a harmonic oscillator of frequency omega_0, constructs the gle kernel for transformation of the velocity velocity autocorrelation function."""
    dw = abs(omega[1] - omega[0])
    ngrid = len(omega)
    dKer = np.zeros((ngrid, ngrid), float)
    omlist = omega.copy()
    omlist[0] = max(omlist[0], dw * 1e-1)
    om2list = omlist ** 2
    y = 0
    if Ap[0, 0] < 2.0 * dw:
        print(
            "# WARNING: White-noise term is weaker than the spacing of the frequency grid. Will increase automatically to avoid instabilities in the numerical integration."
        )

    # outer loop over the physical frequency
    for omega_0 in omlist:
        # works in "scaled coordinates" to stabilize the machinery for small or large omegas
        dAqp = Aqp(omega_0, Ap) / omega_0
        dDqp = Dqp(omega_0, Dp) / omega_0
        dCqp = Cqp(omega_0, dAqp, dDqp)
        if dAqp[1, 1] < 2.0 * dw / omega_0:
            dAqp[1, 1] = 2.0 * dw / omega_0

        dAqp2 = np.dot(dAqp, dAqp)
        # diagonalizes dAqp2 to accelerate the evaluation further down in the inner loop
        w2, omat = np.linalg.eig(dAqp2)
        w = np.sqrt(w2)
        oinv = np.linalg.inv(omat)

        ow1 = omat[1, :] * w / omega_0
        o1cqp1 = np.dot(oinv, dCqp)[:, 1]
        x = 0
        om2om0 = om2list / omega_0 ** 2
        # keeps working in scaled coordinates at this point
        for oo0x in om2om0:
            dKer[x, y] = np.real(np.dot(ow1, o1cqp1 / (w2 + oo0x)))
            x += 1
        y += 1
    return dKer * dw * 2.0 / np.pi


def ISRA(omega, ker, y, dparam, oprefix):
    """Given the thermostatted facf spectrum and the range of frequencies, constructs the vibration density of states"""
    steps = dparam[0]
    stride = dparam[1]
    f = y
    CT = ker.T
    CTC = np.dot(ker.T, ker)
    npad = int(np.log10(steps) + 1)

    for i in range(steps):
        f = f * np.dot(CT, y) / np.dot(CTC, f)
        if np.fmod(i, stride) == 0 and i != 0:
            cnvg = np.asarray(
                (
                    np.linalg.norm((np.dot(f, ker) - y)) ** 2,
                    np.linalg.norm(np.gradient(np.gradient(f))) ** 2,
                )
            )
            dcomm = "# error, laplacian =   " + str(cnvg[0]) + ", " + str(cnvg[1])
            np.savetxt(
                oprefix + "_" + str(i).rjust(npad, "0") + ".data",
                np.vstack((omega, f)).T,
                header=dcomm,
            )
        cnvg = np.asarray(
            (
                i,
                np.linalg.norm((np.dot(f, ker) - y)) ** 2,
                np.linalg.norm(np.gradient(np.gradient(f))) ** 2,
            )
        )
    return f


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
        # opens & parses the i-pi input file
        ifile = open(path2ixml, "r")
        xmlrestart = io_xml.xml_parse_file(ifile)
        ifile.close()

        isimul = InputSimulation()
        isimul.parse(xmlrestart.fields[0][1])
        simul = isimul.fetch()

        # parses the drift and diffusion matrices of the GLE thermostat.
        ttype = str(type(simul.syslist[0].motion.thermostat).__name__)
        P = float(simul.syslist[0].init.nbeads)
        kbT = float(simul.syslist[0].ensemble.temp) * P
        simul.syslist[0].motion.thermostat.temp = kbT

        if ttype == "ThermoGLE":
            Ap = simul.syslist[0].motion.thermostat.A * tscale
            Cp = simul.syslist[0].motion.thermostat.C / kbT
            Dp = np.dot(Ap, Cp) + np.dot(Cp, Ap.T)
        elif ttype == "ThermoNMGLE":
            Ap = simul.syslist[0].motion.thermostat.A[0] * tscale
            Cp = simul.syslist[0].motion.thermostat.C[0] / kbT
            Dp = np.dot(Ap, Cp) + np.dot(Cp, Ap.T)
        elif ttype == "ThermoLangevin" or ttype == "ThermoPILE_L":
            Ap = (
                np.asarray([1.0 / simul.syslist[0].motion.thermostat.tau]).reshape(
                    (1, 1)
                )
                * tscale
            )
            Cp = np.asarray([1.0]).reshape((1, 1))
            Dp = np.dot(Ap, Cp) + np.dot(Cp, Ap.T)
        else:
            raise Exception(
                "GLE thermostat not found. The allowed thermostats are gle, nm_gle, langevin and pile_l."
            )

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
    print("# computing the kernel.")
    ker = gleKernel(ix, Ap, Dp)

    # (de-)convolutes the spectrum
    if action == "conv":
        print("# printing the output spectrum.")
        output_facf((ix, np.dot(iy, ker.T)), oprefix, input_facf(path2ifacf, nrows, 1))
    elif action == "deconv":
        print("# deconvoluting the input spectrum.")
        ISRA(ix, ker, iy, dparam, oprefix)


if __name__ == "__main__":
    # adds description of the program.
    parser = argparse.ArgumentParser(
        description="Given the parameters of a Generalized Langevin Equation and the microcanonical density of states predicts the power spectrum of the velocity-velocity auto-correlation obtained by the dynamics. Conversely, given the power spectrum velocity-velocity autocorrelation function removes the disturbance affected by the thermostat and returns the underlying vibrational density of states. "
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
