#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml


# Disable
def blockPrint():
    sys.stdout = open(os.devnull, "w")


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def get_A(path2iipi):
    """
    Parses the i-PI input to read the relevant data files
    and returns the Helmholtz free energy obtained within the
    self consitent phonons (SCP) approximation.
    """

    blockPrint()

    # Parses the i-PI input file.
    ifile = open(path2iipi, "r")
    xmlrestart = io_xml.xml_parse_file(ifile)
    ifile.close()

    # Initializes the simulation class.
    isimul = InputSimulation()
    isimul.parse(xmlrestart.fields[0][1])
    simul = isimul.fetch()

    enablePrint()

    # Obtains various parameters and relevant filenames.
    prefix = simul.outtemplate.prefix + "." + simul.syslist[0].motion.prefix
    max_iter = simul.syslist[0].motion.max_iter
    batch_exp = simul.syslist[0].motion.batch_weight_exponent
    kbT = float(simul.syslist[0].ensemble.temp)
    beta = 1.0 / kbT

    # Checks for ASR type to skip
    # over zero modes.
    asr = simul.syslist[0].motion.asr
    if asr == "crystal":
        nz = 3
    elif asr == "poly":
        nz = 6

    iD_list = []
    q_list = []
    K_list = []
    V_list = []
    hw_list = []
    x_list = []
    v_list = []

    for i in range(max_iter):
        try:
            # Imports the q, iD, x, f from the i^th  SCP iteration.
            # "0" stands for the SCP trial Hamiltonian.
            # iD -> the covariance matrix of the canonical density in position representation.
            # q  -> the position at minimum potential energy.
            # K  -> the force constant matrix.
            # V  -> the minimum potential energy.
            # hw -> the frequencies.
            # vH -> the ensemble average of the potential w.r.t the canonical density.
            # AH -> the free energy.
            iD_list.append(np.loadtxt(prefix + ".iD." + str(i)))
            q_list.append(np.loadtxt(prefix + ".q." + str(i)))
            K_list.append(np.loadtxt(prefix + ".K." + str(i)))
            V_list.append(np.loadtxt(prefix + ".V0." + str(i)))
            hw_list.append(np.loadtxt(prefix + ".w." + str(i))[nz:])
            # x  -> the samples drawn from the canonical density of the trial distribution.
            # v  -> the potential enegry of the samples.
            x_list.append(np.loadtxt(prefix + ".x." + str(i)))
            v_list.append(np.loadtxt(prefix + ".v." + str(i)))

        except IOError:
            break

    print("# Finished Import")
    for i in range(max_iter):
        try:
            # Imports the q, iD, x, f from the i^th  SCP iteration.
            iD0 = iD_list[i]
            q0 = q_list[i]
            K0 = K_list[i]
            V0 = V_list[i]
            hw0 = hw_list[i]
            betahw0 = beta * hw0
            vH0 = np.sum(hw0 * np.cosh(betahw0 / 2.0) / np.sinh(betahw0 / 2.0) * 0.250)
            AH0 = np.sum(hw0 * 0.5 + kbT * np.log(1 - np.exp(-betahw0)))

        except IOError:
            break

        # Stores the potential /free energy of the initial harmonic trail Hamiltonian.
        if i == 0:
            A_harm = AH0
            v_harm = vH0
            V0_harm = V0
            print(
                (
                    "%23s %23s %23s %23s %23s %23s %23s"
                    % (
                        "# ITERATION",
                        "A_SCP",
                        "A_SCP-A_HARM",
                        "ERROR",
                        "E_SCP",
                        "E_SCP-E_HARM",
                        "ERROR",
                    )
                )
            )

        # Initializes the average and the error in the difference between the physical
        # and the SCP potential.
        avg_dv = 0.0
        err_dv = 0.0
        norm = 0

        # Inner loop over previous SCP steps
        for j in range(i + 1):
            try:
                # Imports the q, iD of the j^th trial Hamiltonian.
                # The idea is to reweight using samples drawn from the j <= i trial Hamiltonians.
                iD = iD_list[j]
                q = q_list[j]

                # x  -> the samples drawn from the canonical density of the trial distribution.
                # v  -> the potential enegry of the samples.
                x = x_list[j]
                v = v_list[j]
            except IOError:
                break

            # vh -> the harmonic component of the potential energy.
            # Note that the harmonic Hamiltonian is the i^th one
            # while the samples are drawn from the j^th one.
            vh = 0.5 * np.sum(np.dot(x - q0, K0.T) * (x - q0), axis=1)

            # Calculates the statistical weight of each sample.
            w = np.exp(
                -(0.50 * np.dot(iD0, (x - q0).T).T * (x - q0)).sum(axis=1)
                + (0.50 * np.dot(iD, (x - q).T).T * (x - q)).sum(axis=1)
            )
            V1 = np.sum(w)
            V2 = np.sum(w**2)

            # Calculates the average amd error (over the samples from the j^th batch of samples)
            # associated with the "anharmonic" component of the potential.
            avg_dv_j = np.nan_to_num(np.sum(w * (v - vh)) / V1)
            err_dv_j = np.nan_to_num(
                np.sum(w * (v - vh - avg_dv_j) ** 2) / (V1 - V2 / V1)
            )

            # Calculates the "batch" weight.
            c = np.nan_to_num(np.exp(-np.var(np.log(w)))) ** batch_exp

            # Accumulates the contribution to the average (and error) from the j^th batch.
            avg_dv += c * avg_dv_j
            err_dv += c**2 * err_dv_j
            norm += c

        avg_dv = np.nan_to_num(avg_dv / norm)
        err_dv = np.nan_to_num(err_dv / norm**2 / len(w))

        # Calculates the SCP potential / free energy.
        A_scp = AH0 + avg_dv
        A_scp_err = np.sqrt(err_dv)
        E_scp = 2.0 * vH0 + avg_dv
        E_scp_err = np.sqrt(err_dv)
        # Calculates the SCP potential / free energy correction
        # w.r.t the initial trial Hamiltonian.
        A_scp_corr = A_scp - A_harm - V0_harm
        E_scp_corr = E_scp - 2.0 * v_harm - V0_harm

        print(
            (
                "%23d %23.8e %23.8e %23.8e %23.8e %23.8e %23.8e"
                % (i, A_scp, A_scp_corr, A_scp_err, E_scp, E_scp_corr, E_scp_err)
            )
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
Calculates the Helmholtz free energy within the self consistent phonons
approximation. The output contains the full self-consistent-phonon
energies, including the minimum energy potential. The center-of-mass
component is not included. The correction relative to the baseline
harmonic description is also reported.
"""
    )
    parser.add_argument(
        "-i",
        "--input_xml",
        type=str,
        required=True,
        help="Path to the i-PI xml file that was used to run the SCP calculation.",
    )
    args = parser.parse_args()
    path2iipi = args.input_xml

    get_A(path2iipi)
