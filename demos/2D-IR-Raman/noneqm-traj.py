#!/usr/bin/env python3

"""Runs nonequilibrium trajectories based on equilibrium-nonequilibrium method for nonlinear spectroscopy."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import os
from math import ceil, floor
import numpy as np
import xml.etree.ElementTree as et
import re

# Check that we have the import path for this i-PI set and if not, add it.
dir_root = os.path.realpath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
)
if not dir_root in sys.path:
    sys.path.insert(0, dir_root)

from ipi.utils.softexit import softexit
from ipi.engine.simulation import Simulation
import ipi.engine.outputs as eoutputs
from ipi.engine.initializer import init_chk
from ipi.utils.messages import verbosity, info


class NonEqmTraj(object):
    """Class containing the details of equilibrium-nonequilibrium 2D spectra simulation.

    Attributes:
       epsilon: Magnitude of the external electric field.
       tsteps: Number of nonequilibrium dynamics steps.
    """

    def __init__(self, epsilon, tsteps):
        """Initializes NonEqmTraj object.

        Args:
           epsilon: Magnitude of the external electric field.
           tsteps: Number of nonequilibrium dynamics steps.
        """
        self.epsilon = epsilon
        self.tsteps = tsteps

    def init_sim_single_step(self, sim):
        sim.tsteps = 1
        outputs = sim.outputs
        sim.outputs = []
        sim.run()
        sim.outputs = outputs

    def fetch_data_and_modify_simulation(self, sim):
        """Gets relevant data from simulation into NonEqmTraj.

        Modifies sim.outputs elements by deleting everything that is not dipole or polarizability. This is
        useful for nonequilibrium trajectories because we might want to print out more data in the equilibrium
        trajectory calculation.

        Also stores useful information like stride and filename of checkpoint/derivative outputs, which we read
        back to initialize nonequilibrium trajectories, as well as the total number of steps of the equilibrium
        trajectory.
        The total number of steps is modified and set to new_tsteps.

        This procedure is invoked once after equilibrium trajectory is computed and before the first nonequilibrium
        trajectory is launched.
        """
        self.nbeads = sim.syslist[0].beads.nbeads
        self.tsteps_eq = sim.tsteps
        self.init_sim_single_step(sim)
        sim.tsteps = self.tsteps
        # Have to loop over an auxiliary list of output elements, because sim.outputs is being modified inside the loop.
        outputs = list(sim.outputs[:])
        for o in outputs:
            if type(o) is eoutputs.TrajectoryOutput:
                if o.what == "extras":
                    if o.extra_type == "polarizability" or o.extra_type == "dipole":
                        continue  # Don't remove this element of output, we want to output dipoles and polarizabilities.
                    if o.extra_type == "dipole_derivative":
                        self.der_fn = self.remove_noneqm_suffix(o.filename)
            # Store values that will help us loop over chk files.
            if type(o) is eoutputs.CheckpointOutput:
                self.chk_stride = o.stride
                self.chk_fn = self.remove_noneqm_suffix(o.filename)
            sim.outputs.remove(
                o
            )  # Remove everything that is not dipole or polarizability.

    def remove_noneqm_suffix(self, string):
        split_string = string.split(".")
        split_string[0] = re.sub("_noneqm$", "", split_string[0])
        return ".".join(split_string)

    def prepare_for_run(self, sim, step, kick):
        """Reads initial q and p from a checkpoint file, applies the kick, and resets step to zero.
        Invoked for each neq trajectory."""
        file_in = self.chk_fn + "_" + str(step)
        new_beads = init_chk(file_in)[0]
        sim.syslist[0].beads.q = new_beads.q
        sim.syslist[0].beads.p = (
            new_beads.p + 0.5 * kick * self.epsilon * self.der[step]
        )
        sim.step = 0

    def run(self, sim, t_first=None, t_last=None):
        """Runs nonequilibrium trajectories."""

        self.fetch_data_and_modify_simulation(sim)

        t_min = ceil(self.tsteps / self.chk_stride)
        t_max = floor(self.tsteps_eq / self.chk_stride)

        t_first = t_min if t_first is None else t_first
        t_last = t_max if t_last is None else t_last
        if (t_first < t_min) or (t_first > t_last):
            sys.exit(
                "t_first should be greater than or equal to "
                + str(t_min)
                + " and less than t_last"
            )
        if t_last > t_max:
            sys.exit("t_last should be less than or equal to " + str(t_max))

        # Bead number formatting with pre-padded zeros (from ipi/engine/outputs.py).
        fmt_bead = (
            "{0:0" + str(int(1 + np.floor(np.log(self.nbeads) / np.log(10)))) + "d}"
        )
        self.der = np.transpose(
            np.array(
                [
                    np.loadtxt(self.der_fn + "_" + fmt_bead.format(b))
                    for b in range(self.nbeads)
                ]
            ),
            [1, 0, 2],
        )
        for step in range(t_first, t_last + 1):
            for kick in [-1, 1]:
                self.prepare_for_run(sim, step, kick)
                sim.run()
                #############################################################################################
                # Stop the thread that monitors for softexit. This is needed because a new thread is started
                # every time we invoke simulation run, leading to a constant increase in the number of
                # threads. The code breaks when a maximum number of threads is reached.
                softexit._doloop[0] = False
                while softexit._thread.is_alive():
                    softexit._thread.join(0.5)
                    # Print out a message saying that this trajectory finished.
                    info(
                        "Trajectory "
                        + str(step)
                        + ["(-)", "(+)"][(kick + 1) // 2]
                        + " finished successfully.",
                        verbosity.quiet,
                    )
                #############################################################################################


def main(fn_input, options):
    """Main procedure:
    1) Modifies fn_input and generates fn_input_noneqm.
    2) Runs nonequilibrium trajectories.
    """

    spec = NonEqmTraj(options.epsilon, options.tsteps)
    tree = et.parse(fn_input)
    root = tree.getroot()
    prefix = root.find("output").attrib["prefix"]
    root.find("output").attrib["prefix"] = str(prefix) + "_noneqm"
    fn_input_noneqm = fn_input + "_noneqm"
    tree.write(fn_input_noneqm)

    simulation = Simulation.load_from_xml(
        open(fn_input_noneqm), request_banner=False, custom_verbosity=options.verbosity
    )
    spec.run(simulation, options.tfirst, options.tlast)
    softexit.trigger(status="success", message=" @ SIMULATION: Exiting cleanly.")


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(
        usage="%prog [options] <input file>",
        description="noneqm-traj runs nonequilibrium trajectories for 2D IR-Raman spectra",
    )

    parser.add_option(
        "-V",
        "--verbosity",
        dest="verbosity",
        default=None,
        choices=["quiet", "low", "medium", "high", "debug"],
        help="Define the verbosity level.",
    )

    parser.add_option(
        "-e",
        "--epsilon",
        dest="epsilon",
        type="float",
        default=0.1,
        help="Epsilon parameter controlling the field strength."
        "For a given epsilon and time step delta_t of dynamics (both in au), the electric field is E = epsilon/delta_t in atomic units of [Eh/(e a0)] = 51.422 [V/Angstrom].",
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

    # Everything is ready. Go!
    main(args[0], options)
