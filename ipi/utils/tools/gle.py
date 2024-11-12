"""
Utility functions to work with generalized Langevin equations
in the "pseudo-Markovian" formulation. See
Ceriotti, Bussi and Parrinello JCTC (2010) http://doi.org/10.1021/ct900563s
for the notation and basic equations
"""

import numpy as np
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml
from ipi.utils.messages import warning, info, verbosity


def get_gle_matrices(ipi_input, time_scaling=1.0):
    """
    Parses an i-PI XML input file and extracts the appropriate GLE
    matrices (Ap, Cp, Dp) from the input. It also does so for Langevin
    dynamics, by taking the tau to determine the friction (app=1/tau)
    """

    # opens & parses the i-pi input file
    ifile = open(ipi_input, "r")
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
        Ap = simul.syslist[0].motion.thermostat.A * time_scaling
        Cp = simul.syslist[0].motion.thermostat.C / kbT
        Dp = np.dot(Ap, Cp) + np.dot(Cp, Ap.T)
    elif ttype == "ThermoNMGLE":
        Ap = simul.syslist[0].motion.thermostat.A[0] * time_scaling
        Cp = simul.syslist[0].motion.thermostat.C[0] / kbT
        Dp = np.dot(Ap, Cp) + np.dot(Cp, Ap.T)
    elif ttype == "ThermoLangevin" or ttype == "ThermoPILE_L":
        Ap = (
            np.asarray([1.0 / simul.syslist[0].motion.thermostat.tau]).reshape((1, 1))
            * time_scaling
        )
        Cp = np.asarray([1.0]).reshape((1, 1))
        Dp = np.dot(Ap, Cp) + np.dot(Cp, Ap.T)
    else:
        raise Exception(
            "GLE thermostat not found. The allowed thermostats are gle, nm_gle, langevin and pile_l."
        )
    return Ap, Cp, Dp


def Aqp(omega_0, Ap):
    """Given the free particle Ap matrix and the frequency
    of a harmonic oscillator, computes the full drift matrix
    for an Ornstein-Uhlenbeck process describing the dynamics
    of (q,p,s)."""

    dAqp = np.zeros(np.asarray(Ap.shape) + 1)
    dAqp[1, 0] = -np.power(omega_0, 2)
    dAqp[0, 1] = 1
    dAqp[1:, 1:] = Ap

    return dAqp


def Dqp(omega_0, Dp):
    """Given the free particle Dp matrix and the frequency of
    the harmonic oscillator, computes the full D matrix
    for an Ornstein-Uhlenbeck process describing the dynamics
    of (q,p,s)."""

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


def Cqp(omega0, dAqp, dDqp):
    """Computes the covariance for (q,p,s), given the drift and squared diffusion matrices
    A and D for the full state vector."""
    # "stabilizes" the calculation by removing the trivial dependence of <a^2> on omega0 until the very end
    dAqp = dAqp.copy()
    dDqp = dDqp.copy()
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


def gle_frequency_kernel(omega, Ap, Dp):
    """Given the Cp and Dp matrices for a harmonic oscillator of frequency omega_0,
    constructs the kernel function that describes the convolution of the velocity
    velocity autocorrelation function of an oscillator of frequency omega_0.
    In other terms, given a list of frequencies omega and free-particle GLE
    matrices Ap and Dp, returns K(w,w') which is the Fourier Transform for
    the v-v correlation of an oscillator of frequency w evaluated at w'."""
    dw = abs(omega[1] - omega[0])
    ngrid = len(omega)
    dKer = np.zeros((ngrid, ngrid), float)
    omlist = omega.copy()
    omlist[0] = max(omlist[0], dw * 1e-1)  # can't handle zero frequency
    om2list = omlist**2
    y = 0
    if Ap[0, 0] < 2.0 * dw:
        warning(
            """White-noise term is weaker than the spacing of the 
            frequency grid. Will increase automatically to avoid 
            instabilities in the numerical integration.""",
            verbosity.low,
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
        om2om0 = om2list / omega_0**2
        # keeps working in scaled coordinates at this point
        for oo0x in om2om0:
            dKer[x, y] = np.real(np.dot(ow1, o1cqp1 / (w2 + oo0x)))
            x += 1
        y += 1

    return dKer * dw * 2.0 / np.pi


def isra_deconvolute(omega, y, ker, steps=100, stride=None):
    """Given the an ACF spectrum computed from a thermostatted trajectories, the range of frequencies, and
    the kernel function that describes the oscillator dynamics for the same GLE parameters,
    reconstructs the vibration density of states of a NVE dynamics of the same system"""

    if stride is None:
        stride = steps // 10
    f = y.copy()
    CT = ker.T
    CTC = np.dot(ker.T, ker)

    spectra = []
    errors = []
    laplas = []
    for i in range(steps):
        f = f * np.dot(CT, y) / np.dot(CTC, f)
        if np.fmod(i, stride) == 0 and i != 0:
            cnvg = np.asarray(
                (
                    np.linalg.norm((np.dot(f, ker) - y)) ** 2,
                    np.linalg.norm(np.gradient(np.gradient(f))) ** 2,
                )
            )
            info(f"# error, laplacian =   {cnvg[0]}, {cnvg[1]}", verbosity.low)
            errors.append(cnvg[0])
            laplas.append(cnvg[1])
            spectra.append(f.copy())

    return f, np.asarray(spectra), np.asarray(errors), np.asarray(laplas)
