"""Creates objects that hold the whole simulation."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from ipi import ipi_global_settings
from ipi.utils.depend import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *
from ipi.utils.prng import *
from ipi.utils.io import *
from ipi.utils.io.inputs.io_xml import *
from ipi.utils.messages import verbosity
from ipi.engine.smotion import Smotion
from ipi.inputs.prng import InputRandom
from ipi.inputs.system import InputSystem, InputSysTemplate
from ipi.engine.system import System
import ipi.inputs.forcefields as iforcefields
import ipi.engine.forcefields as eforcefields
import ipi.inputs.outputs as ioutputs
from ipi.inputs.smotion import InputSmotion


__all__ = ["InputSimulation"]


class InputSimulation(Input):
    """Simulation input class.

    Handles generating the appropriate forcefield class from the xml input file,
    and generating the xml checkpoint tags and data from an instance of the
    object.

    Attributes:
       verbosity: A string saying how much should be output to standard output.
       mode: A string which determines what type of simulation will be run.

    Fields:
       output: A list of the required outputs.
       prng: A random number generator object.
       step: An integer giving the current simulation step. Defaults to 0.
       total_steps: The total number of steps. Defaults to 1000
       total_time:  The wall clock time limit. Defaults to 0 (no limit).

    Dynamic fields:
       system: Holds the data needed to specify the state of a single system.
       ffsocket: Gives a forcefield which will use a socket interface to
          communicate with the driver code.
       fflj: Gives a forcefield which uses the internal Python Lennard-Jones
          script to calculate the potential and forces.
    """

    fields = {
        "prng": (
            InputRandom,
            {
                "help": InputRandom.default_help,
                "default": input_default(factory=Random),
            },
        ),
        "output": (
            ioutputs.InputOutputs,
            {
                "help": ioutputs.InputOutputs.default_help,
                "default": input_default(factory=ioutputs.InputOutputs.make_default),
            },
        ),
        "step": (
            InputValue,
            {"dtype": int, "default": 0, "help": "The current simulation time step."},
        ),
        "total_steps": (
            InputValue,
            {
                "dtype": int,
                "default": 1000,
                "help": "The total number of steps that will be done. If 'step' is equal to or greater than 'total_steps', then the simulation will finish.",
            },
        ),
        "total_time": (
            InputValue,
            {
                "dtype": float,
                "default": 0,
                "help": "The maximum wall clock time (in seconds).",
            },
        ),
        "smotion": (
            InputSmotion,
            {
                "default": input_default(factory=Smotion),
                "help": "Options for a 'super-motion' step between system replicas",
            },
        ),
    }

    attribs = {
        "verbosity": (
            InputAttribute,
            {
                "dtype": str,
                "default": "medium",
                "options": ["quiet", "low", "medium", "high", "debug"],
                "help": "The level of output on stdout.",
            },
        ),
        "threading": (
            InputAttribute,
            {
                "dtype": bool,
                "default": True,
                "help": "Whether multiple-systems execution should be parallel. Makes execution non-reproducible due to the random number generator being used from concurrent threads.",
            },
        ),
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "default": "md",
                "help": "What kind of simulation should be run.",
                "options": ["md", "static"],
            },
        ),
        "safe_stride": (
            InputAttribute,
            {
                "dtype": int,
                "default": 1,
                "help": """Consistent simulation states will be saved every this number of steps. 
Saving state entails a small overhead, so you may want to set this to the smallest output
frequency in your simulation to make i-PI faster. Use at your own risk!
""",
            },
        ),
        "floatformat": (
            InputAttribute,
            {
                "dtype": str,
                "default": "%16.8e",
                "help": "A format for all printed floats.",
            },
        ),
        "sockets_prefix": (
            InputAttribute,
            {
                "dtype": str,
                "default": "/tmp/ipi_",
                "help": "A prefix prepended to the `address` value to form the UNIX-domain socket location.",
            },
        ),
    }

    dynamic = {
        "system": (InputSystem, {"help": InputSystem.default_help}),
        "system_template": (InputSysTemplate, {"help": InputSysTemplate.default_help}),
        "ffsocket": (
            iforcefields.InputFFSocket,
            {"help": iforcefields.InputFFSocket.default_help},
        ),
        "ffdirect": (
            iforcefields.InputFFDirect,
            {"help": iforcefields.InputFFDirect.default_help},
        ),
        "fflj": (
            iforcefields.InputFFLennardJones,
            {"help": iforcefields.InputFFLennardJones.default_help},
        ),
        "ffdmd": (
            iforcefields.InputFFdmd,
            {"help": iforcefields.InputFFdmd.default_help},
        ),
        "ffdebye": (
            iforcefields.InputFFDebye,
            {"help": iforcefields.InputFFDebye.default_help},
        ),
        "ffplumed": (
            iforcefields.InputFFPlumed,
            {"help": iforcefields.InputFFPlumed.default_help},
        ),
        "ffyaff": (
            iforcefields.InputFFYaff,
            {"help": iforcefields.InputFFYaff.default_help},
        ),
        "ffcommittee": (
            iforcefields.InputFFCommittee,
            {"help": iforcefields.InputFFCommittee.default_help},
        ),
        "ffrotations": (
            iforcefields.InputFFRotations,
            {"help": iforcefields.InputFFRotations.default_help},
        ),
        "ffsgdml": (
            iforcefields.InputFFsGDML,
            {"help": iforcefields.InputFFsGDML.default_help},
        ),
        "ffcavphsocket": (
            iforcefields.InputFFCavPhSocket,
            {"help": iforcefields.InputFFCavPhSocket.default_help},
        ),
    }

    default_help = "This is the top level class that deals with the running of the simulation, including holding the simulation specific properties such as the time step and outputting the data."
    default_label = "SIMULATION"

    def store(self, simul):
        """Takes a simulation instance and stores a minimal representation of it.

        Args:
           simul: A simulation object.
        """

        super(InputSimulation, self).store()

        self.sockets_prefix.store(simul.sockets_prefix)
        self.output.store(simul.outtemplate)
        self.prng.store(simul.prng)
        self.step.store(simul.step)
        self.total_steps.store(simul.tsteps)
        self.total_time.store(simul.ttime)
        self.smotion.store(simul.smotion)
        self.threading.store(simul.threading)
        self.floatformat.store(ipi_global_settings["floatformat"])

        # this we pick from the messages class. kind of a "global" but it seems to
        # be the best way to pass around the (global) information on the level of output.
        if verbosity.debug:
            self.verbosity.store("debug")
        elif verbosity.high:
            self.verbosity.store("high")
        elif verbosity.medium:
            self.verbosity.store("medium")
        elif verbosity.low:
            self.verbosity.store("low")
        elif verbosity.quiet:
            self.verbosity.store("quiet")
        else:
            raise ValueError("Invalid verbosity level")

        self.mode.store(simul.mode)
        self.safe_stride.store(simul.safe_stride)

        _fflist = [v for k, v in sorted(simul.fflist.items())]
        if len(self.extra) != len(_fflist) + len(simul.syslist):
            self.extra = [0] * (len(_fflist) + len(simul.syslist))

        for (
            _ii,
            _obj,
        ) in enumerate(_fflist + simul.syslist):
            if self.extra[_ii] == 0:
                if type(_obj) is eforcefields.FFSocket:
                    _iobj = iforcefields.InputFFSocket()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffsocket", _iobj)
                elif isinstance(_obj, eforcefields.FFDirect):
                    _iobj = iforcefields.InputFFDirect()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffdirect", _iobj)
                elif isinstance(_obj, eforcefields.FFLennardJones):
                    _iobj = iforcefields.InputFFLennardJones()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("fflj", _iobj)
                elif isinstance(_obj, eforcefields.FFdmd):
                    _iobj = iforcefields.InputFFdmd()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffdmd", _iobj)
                elif isinstance(_obj, eforcefields.FFDebye):
                    _iobj = iforcefields.InputFFDebye()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffdebye", _iobj)
                elif isinstance(_obj, eforcefields.FFPlumed):
                    _iobj = iforcefields.InputFFPlumed()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffplumed", _iobj)
                elif isinstance(_obj, eforcefields.FFYaff):
                    _iobj = iforcefields.InputFFYaff()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffyaff", _iobj)
                elif isinstance(_obj, eforcefields.FFsGDML):
                    _iobj = iforcefields.InputFFsGDML()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffsgdml", _iobj)
                elif isinstance(_obj, eforcefields.FFCommittee):
                    _iobj = iforcefields.InputFFCommittee()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffcommittee", _iobj)
                elif isinstance(_obj, eforcefields.FFRotations):
                    _iobj = iforcefields.InputFFRotations()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffrotations", _iobj)
                elif isinstance(_obj, eforcefields.FFCavPhSocket):
                    _iobj = iforcefields.InputFFCavPhSocket()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffcavphsocket", _iobj)
                elif isinstance(_obj, System):
                    _iobj = InputSystem()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("system", _iobj)
            else:
                self.extra[_ii][1].store(_obj)

    def fetch(self):
        """Creates a simulation object.

        Returns:
           A simulation object of the appropriate type and with the appropriate
           properties and other objects given the attributes of the
           InputSimulation object.

        Raises:
           TypeError: Raised if one of the file types in the stride keyword
              is incorrect.
        """

        super(InputSimulation, self).fetch()

        # We fetch format and store it in the global variable
        ipi_global_settings["floatformat"] = self.floatformat.fetch()
        try:
            _ = ipi_global_settings["floatformat"] % 1.0
        except:
            print("Error: <simulation> has invalid floatformat attribute.")
            exit(-1)

        # small hack: initialize here the verbosity level -- we really assume to have
        # just one simulation object
        verbosity.level = self.verbosity.fetch()

        syslist = []
        fflist = []
        for k, v in self.extra:
            if k == "system":
                syslist.append(v.fetch())
            elif k == "system_template":
                # This will actually generate automatically a bunch
                # of system objects with the desired properties set
                # automatically to many values.
                syslist += v.fetch()
            elif k in [
                "ffsocket",
                "ffdirect",
                "fflj",
                "ffdebye",
                "ffdmd",
                "ffplumed",
                "ffsgdml",
                "ffyaff",
                "ffcommittee",
                "ffrotations",
                "ffcavphsocket",
            ]:
                new_ff = v.fetch()
                if k == "ffsocket":
                    # overrides ffsocket prefix
                    new_ff.socket.sockets_prefix = self.sockets_prefix.fetch()
                fflist.append(new_ff)

        # this creates a simulation object which gathers all the little bits
        import ipi.engine.simulation as esimulation  # import here as otherwise this is the mother of all circular imports...

        rsim = esimulation.Simulation(
            mode=self.mode.fetch(),
            syslist=syslist,
            fflist=fflist,
            outputs=self.output.fetch(),
            prng=self.prng.fetch(),
            smotion=self.smotion.fetch(),
            step=self.step.fetch(),
            tsteps=self.total_steps.fetch(),
            ttime=self.total_time.fetch(),
            threads=self.threading.fetch(),
            safe_stride=self.safe_stride.fetch(),
            sockets_prefix=self.sockets_prefix.fetch(),
        )

        return rsim
