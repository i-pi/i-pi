#!/usr/bin/env python3

"""Runs equilibrium and nonequilibrium trajectories based on equilibrium-nonequilibrium method for nonlinear spectroscopy.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import os
from copy import deepcopy
from math import ceil, floor
import numpy as np
import xml.etree.ElementTree as et
import re

# Check that we have the import path for this i-PI set and if not, add it.
dir_root = os.path.realpath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
if not dir_root in sys.path:
    sys.path.insert(0, dir_root)

from ipi.utils.softexit import softexit
from ipi.engine.simulation import Simulation
import ipi.engine.outputs as eoutputs
from ipi.engine.initializer import init_chk
from ipi.engine.motion.constrained_dynamics import NVEConstrainedIntegrator
from ipi.utils.messages import verbosity

class NonEqmTraj(object):

    """Class containing the details of equilibrium-nonequilibrium 2D spectra simulation.

    Attributes:
       epsilon: Magnitude of the external electric field.
       tsteps: Number of nonequilibrium dynamics steps.
    """
    @staticmethod
    def load_from_xml(file_in):
        """Loads input file `file_in` and constructs an object of class NonEqmTraj to store the 
           relevant information.
        """ 
        tree = et.parse(file_in)
        root = tree.getroot()
        tsteps = int(root.find('./corr_steps').text)
        epsilon = float(root.find('./epsilon').text)
        return NonEqmTraj(epsilon, tsteps)

    def __init__(self, epsilon=0.1, tsteps=100):
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
           useful for nonequilibrium trajectories because we might want to print out more data from in the equilibrium
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
        #Have to loop over an auxiliary list of output elements, because sim.outputs is being modified inside the loop.
        outputs = list(sim.outputs[:])
        for o in outputs:
            if (type(o) is eoutputs.TrajectoryOutput):
                if o.what == "extras":
                    if (o.extra_type == "polarizability" or o.extra_type == "dipole"):
                        if (o.extra_type == "polarizability"):
                            self.pol_stride = o.stride
                            self.pol_fn = self.remove_noneqm_suffix(o.filename)
                        elif (o.extra_type == "dipole"):
                            self.dip_stride = o.stride
                            self.dip_fn = self.remove_noneqm_suffix(o.filename)
                        continue #Don't remove this element of output, we want to output dipoles and polarizabilities.
                    if o.extra_type == "dipole_derivative":
                        self.der_stride = o.stride
                        self.der_fn = self.remove_noneqm_suffix(o.filename)
            #Store values that will help us loop over chk files.
            if (type(o) is eoutputs.CheckpointOutput):
                self.chk_stride = o.stride
                self.chk_fn =  self.remove_noneqm_suffix(o.filename)
            sim.outputs.remove(o) #Remove everything that is not dipole or polarizability.

    def remove_noneqm_suffix(self, string):
        split_string = string.split('.')
        split_string[0] = re.sub('_noneqm$', '', split_string[0])
        return '.'.join(split_string)

    def prepare_for_run(self, sim, step, kick):
        """Reads initial q and p from a checkpoint file, applies the kick, and resets step to zero. 
           Invoked for each neq trajectory."""
        file_in = self.chk_fn + '_' + str(step)
        new_beads = init_chk(file_in)[0]
        sim.syslist[0].beads.q = new_beads.q
        sim.syslist[0].beads.p = new_beads.p + 0.5 * kick * self.epsilon * self.der[step]
        sim.step = 0

    def run(self, sim):
        """Runs nonequilibrium trajectories."""

        self.fetch_data_and_modify_simulation(sim)
        #Bead number formatting with pre-padded zeros (from ipi/engine/outputs.py).
        fmt_bead = (
            "{0:0"
            + str(int(1 + np.floor(np.log(self.nbeads) / np.log(10))))
            + "d}"
        )
        self.der = np.transpose(np.array([np.loadtxt(self.der_fn + '_' + fmt_bead.format(b)) for b in range(self.nbeads)]), [1,0,2])
        for step in range(ceil(self.tsteps/self.chk_stride), floor(self.tsteps_eq/self.chk_stride) + 1):
            for kick in [-1, 1]:
                self.prepare_for_run(sim, step, kick)
                sim.run()
                #############################################################################################
                #Stop the thread that monitors for softexit. This is needed because a new thread is started
                #every time we invoke simulation run, leading to a constant increase in the number of
                #threads. The code breaks when a maximum number of threads is reached.
                softexit._doloop[0] = False
                while softexit._thread.is_alive():
                    softexit._thread.join(0.5)
                #############################################################################################

def main(fn_input, fn_spec_input, options):
    """Main procedure:
       1) Modifies fn_input and generates fn_input_noneqm.
       2) Runs nonequilibrium trajectories.
       3) Closes the sockets and exits.
    """ 

    spec = NonEqmTraj.load_from_xml(fn_spec_input)
    tree = et.parse(fn_input)
    root = tree.getroot()
    prefix = root.find('output').attrib['prefix']
    root.find('output').attrib['prefix'] = str(prefix) + '_noneqm'
    tree.write(fn_input + '_noneqm')
    
    simulation = Simulation.load_from_xml(fn_input + '_noneqm', request_banner=False, custom_verbosity='quiet')
    spec.run(simulation)
    softexit.trigger(" @ SIMULATION: Exiting cleanly.")


if __name__ == '__main__':

    # TODO: Use argparse once we move to Python 2.7.

    from optparse import OptionParser

    parser = OptionParser(usage='%prog [options] <input file>',
                          description='The main i-PI executable used to run '
                                      'a simulation, given an XML input file.'
                          )

    parser.add_option('-V', '--verbosity', dest='verbosity', default='quiet',
                      choices=['quiet', 'low', 'medium', 'high', 'debug'],
                      help='Define the verbosity level.')

    options, args = parser.parse_args()

    # make sure that we have exactly two input files and that they exist
    if len(args) == 0:
        parser.error('No input file name provided.')
    elif len(args) == 1 or len(args) > 2:
        parser.error('Provide two input file names: one for i-pi, one for spectra.')
    else:
        for fn_in in args:
            if not os.path.exists(fn_in):
                parser.error('Input file not found: {:s}'.format(fn_in))

    # Everything is ready. Go!
    main(args[0], args[1], options)
