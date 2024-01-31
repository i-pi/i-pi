#!/usr/bin/env python3
import os
import numpy as np
import xml.etree.ElementTree as et


def read_xml(file_in):
    tree = et.parse(file_in)
    root = tree.getroot()
    return tree, root


class Initializer:
    """Helper class that serves as a container for input parameters read from input file and command line."""

    @staticmethod
    def load_from_xml(file_in):
        """Loads i-pi input file `file_in`"""
        # ----------------------------
        # I-PI input parameters.
        # ----------------------------
        tree, root = read_xml(file_in)
        # Prefix for the i-pi output files.
        try:
            sim_name = root.find("./output").attrib["prefix"]  # Read simulation name.
        except:
            sim_name = "simulation"
        # Number of beads.
        try:
            nbeads = int(
                root.find("./system/initialize").attrib["nbeads"]
            )  # Read number of beads.
        except:
            nbeads = 1
        # Number of steps between nonequilibrium trajectories (taken from stride of chk files).
        step_traj = int(root.find("./output/checkpoint").attrib["stride"])
        # Stride for printing out dipole/polarizability values.
        for out in root.iter("trajectory"):
            if out.attrib["filename"] in ["dip", "pol"]:
                step_print = int(out.attrib["stride"])
        # Size of dynamics timestep.
        delta_t = float(root.find("./system/motion/dynamics/timestep").text)
        # Ensemble temperature in K converted to beta in atomic units.
        beta = float(root.find("./system/ensemble/temperature").text)
        beta = 1.0 / (3.167e-6 * beta)  # beta in atomic units.

        return Initializer(sim_name, nbeads, step_traj, step_print, beta, delta_t)

    def __init__(self, sim_name, nbeads, step_traj, step_print, beta, delta_t):
        self.sim_name = sim_name
        self.nbeads = nbeads
        self.step_traj = step_traj
        self.step_print = step_print
        self.beta = beta
        self.delta_t = delta_t


class ResponseFunction:
    """Implements two-dimensional equilibrium-nonequilibrium response function.
    The result is (beta / epsilon) <(op2_{+}(t2) - op2_{-}(t2)) op1(-t1)>, where op1 = op[0], op2 = op[1],
    and op is given in input and has two components that are either dip (dipole) or pol (poalrizability).
    Dipole and polarizability are projected on the field polarization given by the field_pol parameter.
    Output is the response function of two-dimensional IR-Raman spectroscopy.
    """

    def __init__(self, input_param):
        self.sim_name = input_param.sim_name
        self.out_name = (
            input_param.out_name + "_" if input_param.out_name is not None else ""
        )
        self.nbeads = input_param.nbeads
        self.beta = input_param.beta
        # Dynamics time step.
        self.delta_t = input_param.delta_t
        # Field strength, defines the magintude of the kick in momentum applied for nonequilibrium trajectories.
        self.epsilon = input_param.epsilon
        # Operators (observables) that will be used to compute the response function.
        self.op = [i.strip() for i in input_param.op.split(",")]
        # Defines which components of dipole/polarizability to use.
        self.field_pol = [int(i) for i in input_param.field_pol.split(",")]
        # Number of nonequilibrium dynamics steps printed out.
        self.tsteps = input_param.tsteps // input_param.step_print
        # Number of equilibrium dynamics steps printed out between two nonequilibrium events.
        self.step_traj = input_param.step_traj // input_param.step_print
        # Stride for printing out dipoles and polarizabilities.
        self.step_print = input_param.step_print
        # First checkpoint file that will be used for computing the response function. Should correspond to a time step greater than tsteps.
        self.chk_first = input_param.tfirst
        # Last checkpoint file that will be used for computing the response function.
        self.chk_last = input_param.tlast

        self.fmt_bead = (
            "{0:0" + str(int(1 + np.floor(np.log(self.nbeads) / np.log(10)))) + "d}"
        )

    def process(self):
        # Get values of dip/pol along field_pol direction and average over beads.
        x = []
        op_index = 0
        for op, field_pol in zip(self.op, self.field_pol):
            x_single = 0
            for bead in range(self.nbeads):
                if op_index == 0:
                    x_single += np.loadtxt(
                        self.sim_name + "." + op + "_" + self.fmt_bead.format(bead)
                    )[:, field_pol]
                else:
                    x_single += np.loadtxt(
                        self.sim_name
                        + "_noneqm."
                        + op
                        + "_"
                        + self.fmt_bead.format(bead)
                    )[:, field_pol]
            op_index += 1
            x_single /= self.nbeads
            x.append(x_single)

        # Compute result:
        c = 0  # Count nonequilibrium events.
        r_t = 0
        i_step = self.tsteps + 1  # Number of steps for t1 and t2.
        first = (
            self.chk_first * self.step_traj
            if self.chk_first is not None
            else self.step_traj
        )
        last = (
            (self.chk_last + 1) * self.step_traj
            if self.chk_last is not None
            else len(x[0])
        )
        for i in range(
            first, last, self.step_traj
        ):  # Start at the first available checkpoint from eq dynamics.
            if i + 1 >= i_step:  # Check if we have enough dynamics to go backward -t1.
                # Backward equilibrium trajectory, first operator.
                eq = np.flip(
                    x[0][i - i_step + 1 : i + 1]
                )  # Equilibrium dynamics part from time t (index i) to t-t1.
                # Difference between forward nonequilibrium trajectories, second operator.
                neq = (
                    x[1][2 * c * i_step : (2 * c + 1) * i_step]
                    - x[1][(2 * c + 1) * i_step : 2 * (c + 1) * i_step]
                )  # Neq dynamics.
                # Compute r_t[j, k] = eq[j] * neq[k].
                r_t += np.outer(eq, neq)
                c += 1
        r_t *= self.beta / (c * self.epsilon)
        r_t = np.gradient(r_t, self.step_print * self.delta_t, axis=0, edge_order=2)

        fn_out = self.out_name + self.op[0] + "_" + self.op[1] + "_neq_2d.dat"
        np.savetxt(
            fn_out,
            r_t,
            fmt="%.5e",
            header="Two-time response function. Time step: "
            + str(self.step_print * self.delta_t)
            + " au. Rows: t1 time. Columns: t2 time. Units: au.",
        )

        print("Response function computed successfully. Generated file " + fn_out + ".")


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(
        usage="%prog [options] <input file>",
        description="noneqm-response comutes response function of "
        "2D IR-Raman spectra from equilibrium and nonequilibrium trajectories. "
        "Rows of the output file correspond to different t1 times, columns correspond to t2 times.",
    )
    parser.add_option(
        "-e",
        "--epsilon",
        dest="epsilon",
        type="float",
        default=0.1,
        help="Epsilon parameter controlling the field strength.",
    )
    parser.add_option(
        "-t",
        "--tsteps",
        dest="tsteps",
        type="int",
        default=100,
        help="Number of time steps for nonequilibrium trajectories.",
    )
    parser.add_option(
        "-b",
        "--op",
        dest="op",
        type="str",
        default="dip, pol",
        help="Two operators (dip or pol) in the form: <op1>, <op2>.",
    )
    parser.add_option(
        "-p",
        "--field_pol",
        dest="field_pol",
        type="str",
        default="2, 8",
        help="Polarization components of operators in op. "
        "For dipole: x=0, y=1, z=2. "
        "For polarizability: xx = 0, xy=1, xz=2, yx=3, ..., zz=8",
    )
    parser.add_option(
        "-i",
        "--tfirst",
        dest="tfirst",
        type="int",
        default=None,
        help="Index of first checkpoint file of equilibrium trajectory that will be used.",
    )
    parser.add_option(
        "-f",
        "--tlast",
        dest="tlast",
        type="int",
        default=None,
        help="Index of last checkpoint file of equilibrium trajectory that will be used.",
    )
    parser.add_option(
        "-o",
        "--out_name",
        dest="out_name",
        type="str",
        default=None,
        help="Prefix for output. "
        "Default is <op1>_<op2>_2d_neq.dat, otherwise <out_name>_<op1>_<op2>_2d_neq.dat",
    )
    options, args = parser.parse_args()

    # make sure that we have exactly two input files and that they exist
    if len(args) == 0:
        parser.error("No input file name provided.")
    elif len(args) > 1:
        parser.error("Provide only one input file name.")
    else:
        for fn_in in args:
            if not os.path.exists(fn_in):
                parser.error("Input file not found: {:s}".format(fn_in))

    # Read parameters from the input xml file of equilibrium dynamics.
    init_param = Initializer.load_from_xml(args[0])
    # Input parametrs from the command line.
    init_param.epsilon = options.epsilon
    init_param.tsteps = options.tsteps
    init_param.op = options.op
    init_param.field_pol = options.field_pol
    init_param.tfirst = options.tfirst
    init_param.tlast = options.tlast
    init_param.out_name = options.out_name

    response = ResponseFunction(init_param)
    response.process()
