"""A class that runs the simulation and outputs results.

The root class for the whole simulation. Contains references to all the top
level objects used in the simulation, and controls all the steps that are
not inherently system dependent, like the running of each time step,
choosing which properties to initialise, and which properties to output.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import tracemalloc
import os
import time
from copy import deepcopy

from ipi.utils.depend import depend_value, dpipe, dproperties
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.messages import verbosity, info, warning, banner
from ipi.utils.softexit import softexit
import ipi.engine.outputs as eoutputs
import ipi.inputs.simulation as isimulation

from concurrent.futures import ThreadPoolExecutor

__all__ = ["Simulation"]


class Simulation:
    """Main simulation object.

    Contains all the references and the main dynamics loop. Also handles the
    initialisation and output.

    Attributes:
        prng: A random number generator object.
        mode: A string specifying what kind of simulation will be run.
        fflist: A list of forcefield objects that can be called to compute energy and forces
        syslist: A list of physical systems
        prng: A random number generator.
        tsteps: The total number of steps.
        ttime: The wall clock time (in seconds).
        outtemplate: A template output object to be used to generate the outputs
            list. This will be used for each of the systems.
        outputs: A list of output objects that should be printed during the run
        paratemp: A helper object for parallel tempering simulations
        chk: A checkpoint object which is kept up-to-date in case of emergency exit
        rollback: If set to true, the state of the simulation at the start
            of the step will be output to a restart file rather than
            the current state of the simulation. This is because we cannot
            restart from half way through a step, only from the beginning of a
            step, so this is necessary for the trajectory to be continuous.
        read_only: If set to true, it creates the simulation object but doesn't  initialize/open the sockets

    Depend objects:
        step: The current simulation step.
    """

    @staticmethod
    def load_from_xml(
        fn_input,
        custom_verbosity=None,
        sockets_prefix=None,
        request_banner=False,
        read_only=False,
    ):
        """Load an XML input file and return a `Simulation` object.

        Arguments:
            fn_input (str): Name of the input file.
            custom_verbosity (str): If not `None`, overrides the verbosity
                specified by the input file.
            request_banner (bool): Whether to print the i-PI banner,
                if verbosity is higher than 'quiet'.
            sockets_prefix (str): Use the specified prefix for all Unix domain sockets.
        """

        # parse the file
        xmlrestart = xml_parse_file(open(fn_input))

        # prepare the simulation input object
        input_simulation = isimulation.InputSimulation()

        # check the input and partition it appropriately
        input_simulation.parse(xmlrestart.fields[0][1])

        # override verbosity if requested
        if custom_verbosity is None:
            # Get from the input file
            custom_verbosity = input_simulation.verbosity.fetch()
        input_simulation.verbosity.value = custom_verbosity

        if sockets_prefix is None:
            # Get from the input file
            sockets_prefix = input_simulation.sockets_prefix.fetch()
        input_simulation.sockets_prefix.value = sockets_prefix

        # print banner if not suppressed and simulation verbose enough
        if request_banner and input_simulation.verbosity.value != "quiet":
            banner()

        # create the simulation object
        simulation = input_simulation.fetch()

        # pipe between the components of the simulation
        simulation.bind(read_only)

        # echo the input file if verbose enough
        if verbosity.low:
            print(" # i-PI loaded input file: ", fn_input)
        elif verbosity.medium:
            print(" --- begin input file content ---")
            ifile = open(fn_input, "r")
            for line in ifile.readlines():
                print(line, end=" ")
            print(" ---  end input file content  ---")
            ifile.close()

        return simulation

    def __init__(
        self,
        mode,
        syslist,
        fflist,
        outputs,
        prng,
        smotion=None,
        step=0,
        tsteps=1000,
        ttime=0,
        threads=False,
        safe_stride=1,
        sockets_prefix="ipi_",
    ):
        """Initialises Simulation class.

        Args:
            mode: What kind of simulation is this
            syslist: A list of system objects
            fflist: A list of forcefield objects
            prng: A random number object.
            smotion: A "super-motion" class specifying what to do with different system replicas
            outputs: A list of output objects.
            step: An optional integer giving the current simulation time step.
                Defaults to 0.
            tsteps: An optional integer giving the total number of steps. Defaults
                to 1000.
            ttime: The simulation running time. Used on restart, to keep a
                cumulative total.
        """

        info(" # Initializing simulation object ", verbosity.low)
        self.prng = prng
        self.mode = mode
        self.threading = threads
        self.safe_stride = safe_stride
        self.sockets_prefix = sockets_prefix

        self.syslist = syslist
        for s in syslist:
            s.prng = self.prng  # bind the system's prng to self prng
            s.init.init_stage1(s)

        # TODO - does this have any meaning now that we introduce the smotion class?
        if self.mode == "md" and len(syslist) > 1:
            warning(
                "Multiple systems will evolve independently in a '"
                + self.mode
                + "' simulation."
            )

        self.fflist = {}
        for f in fflist:
            self.fflist[f.name] = f

        self.outtemplate = outputs

        self._step = depend_value(name="step", value=step)
        self.tsteps = tsteps
        self.ttime = ttime
        self.smotion = smotion

        self.chk = None
        self.rollback = True

    def bind(self, read_only=False):
        """Calls the bind routines for all the objects in the simulation."""

        if self.tsteps <= self.step:
            raise ValueError(
                "Simulation has already run for total_steps, will not even start. "
                "Modify total_steps or step counter to continue."
            )

        # initializes the output maker so it can be passed around to systems
        f_start = self.step == 0  # special operations if we're starting from scratch
        if f_start:
            mode = "w"
        else:
            mode = "a"
        self.output_maker = eoutputs.OutputMaker(self.outtemplate.prefix, f_start)

        for s in self.syslist:
            # binds important computation engines
            s.bind(self)

        if read_only:  # returns before we open the sockets
            return

        # start forcefields here so we avoid having a shitload of files printed
        # out only to find the socket is busy or whatever prevented starting the threads
        for k, f in self.fflist.items():
            f.bind(self.output_maker)
            f.start()

        # Checks for repeated filenames.
        filename_list = [x.filename for x in self.outtemplate]
        if len(filename_list) > len(set(filename_list)):
            raise ValueError(
                "Output filenames are not unique. Modify filename attributes."
            )

        self.outputs = []
        for o in self.outtemplate:
            dco = deepcopy(o)  # avoids overwriting the actual filename
            if self.outtemplate.prefix != "":
                dco.filename = self.outtemplate.prefix + "." + o.filename
            if (
                type(dco) is eoutputs.CheckpointOutput
            ):  # checkpoints are output per simulation
                dco.bind(self)
                dpipe(
                    dco._step, o._step
                )  # makes sure that the checkpoint step is updated also in the template
                self.outputs.append(dco)
            else:  # properties and trajectories are output per system
                isys = 0
                for s in self.syslist:  # create multiple copies
                    no = deepcopy(dco)
                    if s.prefix != "":
                        no.filename = s.prefix + "_" + no.filename
                    no.bind(s, mode)
                    self.outputs.append(no)
                    if f_start:  # starting of simulation, print headers (if any)
                        no.print_header()
                    isys += 1

        if self.threading:
            self.executor = ThreadPoolExecutor(
                max_workers=max(len(self.syslist), len(self.outputs))
            )

        self.chk = eoutputs.CheckpointOutput("RESTART", 1, True, 0)
        self.chk.bind(self)

        if self.smotion is not None:
            self.smotion.bind(self.syslist, self.prng, self.output_maker)

    def softexit(self):
        """Deals with a soft exit request.

        Tries to ensure that a consistent restart checkpoint is
        written out.
        """

        if self.step < self.tsteps:
            self.step += 1
        if not self.rollback:
            info("SOFTEXIT: Saving the latest status at the end of the step")
            self.chk.store()

        self.chk.write(store=False)

    def run(self):
        """Runs the simulation.

        Does all the simulation steps, and outputs data to the appropriate files
        when necessary. Also deals with starting and cleaning up the threads used
        in the communication between the driver and the PIMD code.
        """

        # registers the softexit routine
        softexit.register_function(self.softexit)
        softexit.start(self.ttime)

        # starts tracemalloc to debug memory leaks
        if verbosity.debug:
            tracemalloc.start(10)

        # prints inital configuration -- only if we are not restarting
        if self.step == 0:
            self.step = -1
            # must use multi-threading to avoid blocking in multi-system runs with WTE
            if self.threading:
                stepthreads = []
                for o in self.outputs:
                    st = self.executor.submit(o.write)
                    stepthreads.append(st)

                for st in stepthreads:
                    st.result()
            else:
                for o in self.outputs:
                    o.write()
            self.step = 0

        steptime = 0.0
        simtime = time.time()

        cstep = 0
        # tptime = 0.0
        # tqtime = 0.0
        # tttime = 0.0
        ttot = 0.0
        # main MD loop
        for self.step in range(self.step, self.tsteps):
            # stores the state before doing a step.
            # this is a bit time-consuming but makes sure that we can honor soft
            # exit requests without screwing the trajectory

            steptime = -time.time()
            if softexit.triggered:
                break

            # save a consistent state of the simulation that will be saved as a RESTART file in case of premature (soft) exit
            if self.step % self.safe_stride == 0:
                self.chk.store()

            if len(self.syslist) > 0 and self.threading:
                stepthreads = []
                # steps through all the systems
                for s in self.syslist:
                    # creates separate threads for the different systems
                    st = self.executor.submit(s.motion.step, step=self.step)
                    stepthreads.append(st)

                for st in stepthreads:
                    st.result()
            else:
                for s in self.syslist:
                    s.motion.step(step=self.step)

            if softexit.triggered:
                # Don't continue if we are about to exit.
                break

            # does the "super motion" step
            if self.smotion is not None:
                self.smotion.step(self.step)

            if softexit.triggered:
                # Don't write if we are about to exit.
                break

            if self.threading:
                stepthreads = []
                for o in self.outputs:
                    if o.active():  # don't start a thread if it's not needed
                        st = self.executor.submit(o.write)
                        stepthreads.append(st)

                for st in stepthreads:
                    st.result()
            else:
                for o in self.outputs:
                    o.write()

            steptime += time.time()
            ttot += steptime
            cstep += 1

            if (
                verbosity.high
                or (verbosity.medium and self.step % 100 == 0)
                or (verbosity.low and self.step % 1000 == 0)
            ):
                info(
                    " # Average timings at MD step % 7d. t/step: %10.5e"
                    % (self.step, ttot / cstep)
                )
                cstep = 0
                ttot = 0.0

                # tracemalloc memory traces
                if verbosity.debug:
                    snapshot = tracemalloc.take_snapshot()
                    top_stats = snapshot.statistics("lineno")
                    info(" # DEBUG # Top 10 memory allocations: ")
                    for stat in top_stats[:10]:
                        info(stat)
                    top_stats = snapshot.statistics("traceback")
                    info(" # DEBUG # Trace of the top memory allocation:")
                    for line in top_stats[0].traceback.format():
                        info(line)

            if os.path.exists("EXIT"):
                info(" # EXIT file detected! Bye bye!", verbosity.low)
                break

            if (self.ttime > 0) and (time.time() - simtime > self.ttime):
                info(" # Wall clock time expired! Bye bye!", verbosity.low)
                break

        self.rollback = False


dproperties(Simulation, ["step"])
