"""Creates objects that deal with the evaluation of interactions."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from copy import copy
import numpy as np

from ipi.engine.forcefields import (
    ForceField,
    FFSocket,
    FFDirect,
    FFLennardJones,
    FFDebye,
    FFPlumed,
    FFYaff,
    FFsGDML,
    FFCommittee,
    FFdmd,
    FFCavPhSocket,
    FFRotations,
)
from ipi.interfaces.sockets import InterfaceSocket
from ipi.pes import __drivers__
import ipi.engine.initializer
from ipi.inputs.initializer import *
from ipi.utils.inputvalue import *
from ipi.utils.messages import verbosity, warning
from ipi.utils.prng import Random
from ipi.inputs.prng import InputRandom

__all__ = [
    "InputFFSocket",
    "InputFFDirect",
    "InputFFLennardJones",
    "InputFFDebye",
    "InputFFPlumed",
    "InputFFYaff",
    "InputFFsGDML",
    "InputFFCommittee",
    "InputFFdmd",
    "InputFFCavPhSocket",
    "InputFFRotations",
]


class InputForceField(Input):
    """ForceField input class.

    Handles generating one instance of a particular forcefield class from the xml
    input file, and generating the xml checkpoint tags and data from an
    instance of the object.

    Attributes:
       name: The name by which the forcefield will be identified in the System forces section.
       pbc: A boolean describing whether periodic boundary conditions will
          be applied to the atom positions before they are sent to the driver
          code. Default is False.

    Fields:
       latency: The number of seconds to sleep between looping over the requests.
       parameters: A dictionary containing the forcefield parameters.
       activelist: A list of indexes (starting at 0) of the atoms that will be active in this force field.
    """

    attribs = {
        "name": (
            InputAttribute,
            {
                "dtype": str,
                "help": "Mandatory. The name by which the forcefield will be identified in the System forces section.",
            },
        ),
        "pbc": (
            InputAttribute,
            {
                "dtype": bool,
                "default": False,
                "help": "Applies periodic boundary conditions to the atoms coordinates before passing them on to the driver code.",
            },
        ),
        "threaded": (
            InputAttribute,
            {
                "dtype": bool,
                "default": False,
                "help": "Whether the forcefield should use a thread loop to evaluate, or work in serial",
            },
        ),
    }
    fields = {
        "latency": (
            InputValue,
            {
                "dtype": float,
                "default": 1e-4,
                "help": "The number of seconds the polling thread will wait between exhamining the list of requests.",
            },
        ),
        "offset": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "dimension": "energy",
                "help": "A constant offset that is subtracted from the forcefield energy. Useful when there is a large core energy contribution that is constant throughout a simulation and hides significant changes in the 10th digit.",
            },
        ),
        "parameters": (
            InputValue,
            {"dtype": dict, "default": {}, "help": "The parameters of the force field"},
        ),
        "activelist": (
            InputArray,
            {
                "dtype": int,
                "default": np.array([-1]),
                #                                     "default" : input_default(factory=np.array, args =[-1]),
                "help": "List with indexes of the atoms that this socket is taking care of.    Default: [-1] (corresponding to all)",
            },
        ),
    }

    default_help = "Base forcefield class that deals with the assigning of force calculation jobs and collecting the data."
    default_label = "FORCEFIELD"

    def store(self, ff):
        """Takes a ForceField instance and stores a minimal representation of it.

        Args:
           forceb: A ForceField object.
        """

        Input.store(self, ff)
        self.name.store(ff.name)
        self.latency.store(ff.latency)
        self.offset.store(ff.offset)
        self.parameters.store(ff.pars)
        self.pbc.store(ff.dopbc)
        self.activelist.store(ff.active)
        self.threaded.store(ff.threaded)

    _FFCLASS = ForceField

    def fetch(self):
        """Creates a ForceField object.

        Returns:
           A ForceField object.
        """

        super(InputForceField, self).fetch()

        return self._FFCLASS(
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            active=self.activelist.fetch(),
            threaded=self.threaded.fetch(),
        )


class InputFFSocket(InputForceField):
    """Creates a ForceField object with a socket interface.

    Handles generating one instance of a socket interface forcefield class.

    Attributes:
       mode: Describes whether the socket will be a unix or an internet socket.

    Fields:
       address: The server socket binding address.
       port: The port number for the socket.
       slots: The number of clients that can queue for connections at any one
          time.
       timeout: The number of seconds that the socket will wait before assuming
          that the client code has died. If 0 there is no timeout.
    """

    fields = {
        "address": (
            InputValue,
            {
                "dtype": str,
                "default": "localhost",
                "help": "This gives the server address that the socket will run on.",
            },
        ),
        "port": (
            InputValue,
            {
                "dtype": int,
                "default": 65535,
                "help": "This gives the port number that defines the socket.",
            },
        ),
        "slots": (
            InputValue,
            {
                "dtype": int,
                "default": 4,
                "help": "This gives the number of client codes that can queue at any one time.",
            },
        ),
        "exit_on_disconnect": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "Determines if i-PI should quit when a client disconnects.",
            },
        ),
        "timeout": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "This gives the number of seconds before assuming a calculation has died. If 0 there is no timeout.",
            },
        ),
    }
    attribs = {
        "mode": (
            InputAttribute,
            {
                "dtype": str,
                "options": ["unix", "inet"],
                "default": "inet",
                "help": "Specifies whether the driver interface will listen onto a internet socket [inet] or onto a unix socket [unix].",
            },
        ),
        "matching": (
            InputAttribute,
            {
                "dtype": str,
                "options": ["auto", "any", "lock"],
                "default": "auto",
                "help": "Specifies whether requests should be dispatched to any client, automatically matched to the same client when possible [auto] or strictly forced to match with the same client [lock].",
            },
        ),
    }

    attribs.update(InputForceField.attribs)
    fields.update(InputForceField.fields)

    # FFSocket polling mechanism won't work with non-threaded execution
    attribs["threaded"] = (
        InputValue,
        {
            "dtype": bool,
            "default": True,
            "help": "Whether the forcefield should use a thread loop to evaluate, or work in serial. Should be set to True for FFSockets",
        },
    )

    default_help = "Deals with the assigning of force calculation jobs to different driver codes, and collecting the data, using a socket for the data communication."
    default_label = "FFSOCKET"

    def store(self, ff):
        """Takes a ForceField instance and stores a minimal representation of it.

        Args:
           ff: A ForceField object with a FFSocket forcemodel object.
        """

        if type(ff) not in [FFSocket, FFCavPhSocket]:
            raise TypeError(
                "The type " + type(ff).__name__ + " is not a valid socket forcefield"
            )

        super(InputFFSocket, self).store(ff)

        self.address.store(ff.socket.address)
        self.port.store(ff.socket.port)
        self.timeout.store(ff.socket.timeout)
        self.slots.store(ff.socket.slots)
        self.mode.store(ff.socket.mode)
        self.matching.store(ff.socket.match_mode)
        self.exit_on_disconnect.store(ff.socket.exit_on_disconnect)
        self.threaded.store(True)  # hard-coded

    def fetch(self):
        """Creates a ForceSocket object.

        Returns:
           A ForceSocket object with the correct socket parameters.
        """

        if self.threaded.fetch() is False:
            raise ValueError("FFSockets cannot poll without threaded mode.")
        # just use threaded throughout

        # if using forced match mode, ensure softexit called upon disconnection of a client.
        if self.matching.fetch() == "lock":
            warning(
                'When using matching="lock" pay attention to the possibility of superfluous drivers idling if there are more client codes connected than there are replicas.',
                verbosity.low,
            )
            self.exit_on_disconnect.store(True)

        return FFSocket(
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            active=self.activelist.fetch(),
            threaded=self.threaded.fetch(),
            interface=InterfaceSocket(
                address=self.address.fetch(),
                port=self.port.fetch(),
                slots=self.slots.fetch(),
                mode=self.mode.fetch(),
                timeout=self.timeout.fetch(),
                match_mode=self.matching.fetch(),
                exit_on_disconnect=self.exit_on_disconnect.fetch(),
            ),
        )

    def check(self):
        """Deals with optional parameters."""

        super(InputFFSocket, self).check()
        if self.port.fetch() < 1 or self.port.fetch() > 65535:
            raise ValueError(
                "Port number " + str(self.port.fetch()) + " out of acceptable range."
            )
        elif self.port.fetch() < 1025:
            warning(
                "Low port number being used, this may interrupt important system processes.",
                verbosity.low,
            )

        if self.slots.fetch() < 1 or self.slots.fetch() > 5:
            raise ValueError(
                "Slot number " + str(self.slots.fetch()) + " out of acceptable range."
            )
        if self.latency.fetch() < 0:
            raise ValueError("Negative latency parameter specified.")
        if self.timeout.fetch() < 0.0:
            raise ValueError("Negative timeout parameter specified.")


class InputFFDirect(InputForceField):
    fields = {
        "pes": (
            InputValue,
            {
                "dtype": str,
                "default": "dummy",
                "options": list(__drivers__.keys()),
                "help": "Type of PES that should be used to evaluate the forcefield",
            },
        ),
    }
    fields.update(InputForceField.fields)

    attribs = {}
    attribs.update(InputForceField.attribs)

    default_help = """ Direct potential that evaluates forces through a Python
    call, using PES providers from a list of possible external codes. The available
    PES interfaces are listed into the `ipi/pes` folder, and are the same available
    for the Python driver. The `<parameters>` field should contain a dictionary
    of the specific option of the chosen PES.
    """
    default_label = "FFDirect"

    def store(self, ff):
        super().store(ff)
        self.pes.store(ff.pes)

    def fetch(self):
        super().fetch()

        return FFDirect(
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            threaded=self.threaded.fetch(),
            pes=self.pes.fetch(),
        )


class InputFFLennardJones(InputForceField):
    attribs = {}
    attribs.update(InputForceField.attribs)

    default_help = """Simple, internal LJ evaluator without cutoff, neighbour lists or minimal image convention.
                   Expects standard LJ parameters, e.g. { eps: 0.1, sigma: 1.0 }. """
    default_label = "FFLJ"

    _FFCLASS = FFLennardJones


class InputFFdmd(InputForceField):
    fields = {
        "dmd_coupling": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Specifies the coupling between atom pairs (should be size N*(N-1)/2 ordered c21, c32, c31, c43, c42, c41 etc.  -- in atomic units!)",
            },
        ),
        "dmd_freq": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "Frequency of the oscillation of the time-dependent term",
                "dimension": "frequency",
            },
        ),
        "dmd_dt": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "Time step of the oscillating potential. Should match time step of simulation",
                "dimension": "time",
            },
        ),
        "dmd_step": (
            InputValue,
            {"dtype": int, "default": 0, "help": "The current step counter for dmd."},
        ),
    }

    fields.update(InputForceField.fields)

    attribs = {}
    attribs.update(InputForceField.attribs)

    default_help = """Simple, internal DMD evaluator without without neighbor lists, but with PBC.
                   Expects coupling elements (n*(n-1)/2 of them), oscillating frequency and time step. """
    default_label = "FFDMD"

    def store(self, ff):
        super(InputFFdmd, self).store(ff)
        self.dmd_coupling.store(ff.coupling)
        self.dmd_freq.store(ff.freq)
        self.dmd_dt.store(ff.dtdmd)
        self.dmd_step.store(ff.dmdstep)

    def fetch(self):
        super(InputFFdmd, self).fetch()

        return FFdmd(
            coupling=self.dmd_coupling.fetch(),
            freq=self.dmd_freq.fetch(),
            dtdmd=self.dmd_dt.fetch(),
            dmdstep=self.dmd_step.fetch(),
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            threaded=self.threaded.fetch(),
        )


class InputFFDebye(InputForceField):
    fields = {
        "hessian": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Specifies the Hessian of the harmonic potential. "
                "Default units are atomic. Units can be specified only by xml attribute. "
                r"Implemented options are: 'atomic_unit', 'ev/ang^2'",
                "dimension": "hessian",
            },
        ),
        "x_reference": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Minimum-energy configuration for the harmonic potential",
                "dimension": "length",
            },
        ),
        "v_reference": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "Zero-value of energy for the harmonic potential",
                "dimension": "energy",
            },
        ),
    }

    fields.update(InputForceField.fields)

    attribs = {}
    attribs.update(InputForceField.attribs)

    default_help = """Harmonic energy calculator """
    default_label = "FFDEBYE"

    def store(self, ff):
        super(InputFFDebye, self).store(ff)
        self.hessian.store(ff.H)
        self.x_reference.store(ff.xref)
        self.v_reference.store(ff.vref)

    def fetch(self):
        super(InputFFDebye, self).fetch()

        return FFDebye(
            H=self.hessian.fetch(),
            xref=self.x_reference.fetch(),
            vref=self.v_reference.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            threaded=self.threaded.fetch(),
        )


class InputFFPlumed(InputForceField):
    fields = {
        "file": (
            InputInitFile,
            {
                "default": input_default(
                    factory=ipi.engine.initializer.InitFile, kwargs={"mode": "xyz"}
                ),
                "help": "This describes the location to read the reference structure file from.",
            },
        ),
        "plumed_dat": (
            InputValue,
            {"dtype": str, "default": "plumed.dat", "help": "The PLUMED input file"},
        ),
        "plumed_step": (
            InputValue,
            {
                "dtype": int,
                "default": 0,
                "help": "The current step counter for PLUMED calls",
            },
        ),
        "compute_work": (
            InputValue,
            {
                "dtype": bool,
                "default": True,
                "help": """Compute the work done by the metadynamics bias 
                (to correct the conserved quantity). Note that this might 
                require multiple evaluations of the conserved quantities,
                and can add some overhead.""",
            },
        ),
        "plumed_extras": (
            InputArray,
            {
                "dtype": str,
                "default": input_default(factory=np.zeros, args=(0, str)),
                "help": """List of variables defined in the PLUMED input, 
                that should be transferred to i-PI as `extras` fields. 
                Note that a version of PLUMED greater or equal than 2.10 is necessary 
                to retrieve variables into i-PI, 
                and that, if you are using ring-polymer contraction, the associated force block 
                should use `interpolate_extras` to make sure all beads have values. 
                """,
            },
        ),
    }

    attribs = {}

    attribs.update(InputForceField.attribs)
    fields.update(InputForceField.fields)

    default_help = """ 
Direct PLUMED interface. Can be used to implement metadynamics in i-PI in combination with
the <metad> SMotion class. NB: if you use PLUMED for constant biasing (e.g. for umbrella sampling)
the bias will be computed but there will be no output as specified in the `plumed.dat` file 
unless you include a <metad> tag, that triggers the log update. """

    default_label = "FFPLUMED"

    def store(self, ff):
        super(InputFFPlumed, self).store(ff)
        self.plumed_dat.store(ff.plumed_dat)
        self.plumed_step.store(ff.plumed_step)
        self.compute_work.store(ff.compute_work)
        self.file.store(ff.init_file)
        self.plumed_extras.store(np.array(list(ff.plumed_data.keys())))

    def fetch(self):
        super(InputFFPlumed, self).fetch()

        return FFPlumed(
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            threaded=self.threaded.fetch(),
            compute_work=self.compute_work.fetch(),
            plumed_dat=self.plumed_dat.fetch(),
            plumed_step=self.plumed_step.fetch(),
            plumed_extras=self.plumed_extras.fetch(),
            init_file=self.file.fetch(),
        )


class InputFFYaff(InputForceField):
    fields = {
        "yaffpara": (
            InputValue,
            {
                "dtype": str,
                "default": "parameters.txt",
                "help": "This gives the file name of the Yaff input parameter file.",
            },
        ),
        "yaffsys": (
            InputValue,
            {
                "dtype": str,
                "default": "system.chk",
                "help": "This gives the file name of the Yaff input system file.",
            },
        ),
        "yafflog": (
            InputValue,
            {
                "dtype": str,
                "default": "yaff.log",
                "help": "This gives the file name of the Yaff output log file.",
            },
        ),
        "rcut": (
            InputValue,
            {
                "dtype": float,
                "default": 18.89726133921252,
                "help": "This gives the real space cutoff used by all pair potentials in atomic units.",
            },
        ),
        "alpha_scale": (
            InputValue,
            {
                "dtype": float,
                "default": 3.5,
                "help": "This gives the alpha parameter in the Ewald summation based on the real-space cutoff: alpha = alpha_scale / rcut. Higher values for this parameter imply a faster convergence of the reciprocal terms, but a slower convergence in real-space.",
            },
        ),
        "gcut_scale": (
            InputValue,
            {
                "dtype": float,
                "default": 1.1,
                "help": "This gives the reciprocale space cutoff based on the alpha parameter: gcut = gcut_scale * alpha. Higher values for this parameter imply a better convergence in the reciprocal space.",
            },
        ),
        "skin": (
            InputValue,
            {
                "dtype": int,
                "default": 0,
                "help": "This gives the skin parameter for the neighborlist.",
            },
        ),
        "smooth_ei": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "This gives the flag for smooth truncations for the electrostatic interactions.",
            },
        ),
        "reci_ei": (
            InputValue,
            {
                "dtype": str,
                "default": "ewald",
                "help": "This gives the method to be used for the reciprocal contribution to the electrostatic interactions in the case of periodic systems. This must be one of 'ignore' or 'ewald'. The 'ewald' option is only supported for 3D periodic systems.",
            },
        ),
    }

    fields.update(InputForceField.fields)

    attribs = {}
    attribs.update(InputForceField.attribs)

    default_help = """Uses a Yaff force field to compute the forces."""
    default_label = "FFYAFF"

    def store(self, ff):
        super(InputFFYaff, self).store(ff)
        self.yaffpara.store(ff.yaffpara)
        self.yaffsys.store(ff.yaffsys)
        self.yafflog.store(ff.yafflog)
        self.rcut.store(ff.rcut)
        self.alpha_scale.store(ff.alpha_scale)
        self.gcut_scale.store(ff.gcut_scale)
        self.skin.store(ff.skin)
        self.smooth_ei.store(ff.smooth_ei)
        self.reci_ei.store(ff.reci_ei)

    def fetch(self):
        super(InputFFYaff, self).fetch()

        return FFYaff(
            yaffpara=self.yaffpara.fetch(),
            yaffsys=self.yaffsys.fetch(),
            yafflog=self.yafflog.fetch(),
            rcut=self.rcut.fetch(),
            alpha_scale=self.alpha_scale.fetch(),
            gcut_scale=self.gcut_scale.fetch(),
            skin=self.skin.fetch(),
            smooth_ei=self.smooth_ei.fetch(),
            reci_ei=self.reci_ei.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            threaded=self.threaded.fetch(),
        )


class InputFFsGDML(InputForceField):
    fields = {
        "sGDML_model": (
            InputValue,
            {
                "dtype": str,
                "default": None,
                "help": "This gives the file name of the sGDML model.",
            },
        ),
    }

    fields.update(InputForceField.fields)

    attribs = {}
    attribs.update(InputForceField.attribs)

    default_help = """A SGDML energy calculator """
    default_label = "FFsGDML"

    def store(self, ff):
        super(InputFFsGDML, self).store(ff)
        self.sGDML_model.store(ff.sGDML_model)

    def fetch(self):
        super(InputFFsGDML, self).fetch()

        return FFsGDML(
            sGDML_model=self.sGDML_model.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            threaded=self.threaded.fetch(),
        )


class InputFFCommittee(InputForceField):
    default_help = """Combines multiple forcefields to build a committee model, that can 
                      be used to compute uncertainty-quantified machine-learning models. 
                      Each forcefield can be any of the other FF objects, and each should
                      be used with a client that generates a slightly different estimation
                      of energy and forces. These are averaged, and the mean used as the 
                      actual forcefield. Statistics about the distribution are also returned
                      as extras fields, and can be printed for further postprocessing. 
                      It is also possible for a single FF object to return a JSON-formatted
                      string containing entries `committee_pot`, `committee_force` and 
                      `committee_virial`, that contain multiple members at once. These
                      will be unpacked and combined with whatever else is present. 

                      Also contains options to use it for uncertainty estimation and for
                      active learning in a ML context, based on a committee model.
                      Implements the approaches discussed in 
                      [Musil et al.](http://doi.org/10.1021/acs.jctc.8b00959)
                      and [Imbalzano et al.](http://doi.org/10.1063/5.0036522)
                      """
    default_label = "FFCOMMITTEE"

    dynamic = {
        "ffsocket": (InputFFSocket, {"help": InputFFSocket.default_help}),
        "ffdirect": (InputFFDirect, {"help": InputFFDirect.default_help}),
        "fflj": (InputFFLennardJones, {"help": InputFFLennardJones.default_help}),
        "ffdebye": (InputFFDebye, {"help": InputFFDebye.default_help}),
        "ffplumed": (InputFFPlumed, {"help": InputFFPlumed.default_help}),
        "ffyaff": (InputFFYaff, {"help": InputFFYaff.default_help}),
        "ffsgdml": (InputFFsGDML, {"help": InputFFsGDML.default_help}),
    }

    fields = copy(InputForceField.fields)

    fields["weights"] = (
        InputArray,
        {
            "dtype": float,
            "default": np.array([]),
            "help": """List of weights to be given to the forcefields. Defaults to 1 for each FF. 
                Note that the components are divided by the number of FF, and so the default corresponds to an average.""",
        },
    )
    fields["alpha"] = (
        InputValue,
        {
            "dtype": float,
            "default": 1.0,
            "help": "Scaling of the variance of the model, corresponding to a calibration of the error ",
        },
    )

    fields["baseline_name"] = (
        InputValue,
        {
            "dtype": str,
            "default": "",
            "help": """Name of the forcefield object that should be treated as the baseline for a weighted baseline model.""",
        },
    )

    fields["baseline_uncertainty"] = (
        InputValue,
        {
            "dtype": float,
            "default": -1.0,
            "dimension": "energy",
            "help": """Corresponds to the expected error of the baseline model. This represents the error on the TOTAL potential energy of the simulation. """,
        },
    )
    fields["active_thresh"] = (
        InputValue,
        {
            "dtype": float,
            "default": 0.0,
            "dimension": "energy",
            "help": """The uncertainty threshold for active learning. Structure with an uncertainty above this 
                        value are printed in the specified output file so they can be used for active learning.""",
        },
    )
    fields["active_output"] = (
        InputValue,
        {
            "dtype": str,
            "default": "active_output",
            "help": "Output filename for structures that exceed the accuracy threshold of the model, to be used in active learning.",
        },
    )
    fields["parse_json"] = (
        InputValue,
        {
            "dtype": bool,
            "default": False,
            "help": "Tell the model whether to parse extras string looking for committee values of potential, forces, and virials. Default: false.",
        },
    )

    def store(self, ff):
        """Store all the sub-forcefields"""

        super(InputFFCommittee, self).store(ff)
        _fflist = ff.fflist
        if len(self.extra) != len(_fflist):
            self.extra = [0] * len(_fflist)
        self.weights.store(ff.ffweights)
        self.alpha.store(ff.alpha)
        self.baseline_uncertainty.store(ff.baseline_uncertainty)
        self.baseline_name.store(ff.baseline_name)
        self.active_thresh.store(ff.active_thresh)
        self.active_output.store(ff.active_out)
        self.parse_json.store(ff.parse_json)

        for _ii, _obj in enumerate(_fflist):
            if self.extra[_ii] == 0:
                if isinstance(_obj, FFSocket):
                    _iobj = InputFFSocket()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffsocket", _iobj)
                elif isinstance(_obj, FFLennardJones):
                    _iobj = InputFFLennardJones()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("fflj", _iobj)
                elif isinstance(_obj, FFQUIP):
                    _iobj = InputFFQUIP()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffquip", _iobj)
                elif isinstance(_obj, FFDebye):
                    _iobj = InputFFDebye()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffdebye", _iobj)
                elif isinstance(_obj, FFPlumed):
                    _iobj = InputFFPlumed()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffplumed", _iobj)
                elif isinstance(_obj, FFYaff):
                    _iobj = InputFFYaff()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffyaff", _iobj)
                elif isinstance(_obj, FFsGDML):
                    _iobj = InputFFsGDML()
                    _iobj.store(_obj)
                    self.extra[_ii] = ("ffsgdml", _iobj)
            else:
                self.extra[_ii][1].store(_obj)

    def fetch(self):
        """Fetches all of the FF objects"""

        super(InputFFCommittee, self).fetch()

        fflist = []
        for k, v in self.extra:
            fflist.append(v.fetch())

        # TODO: will actually need to create a FF object here!
        return FFCommittee(
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            fflist=fflist,
            ffweights=self.weights.fetch(),
            alpha=self.alpha.fetch(),
            baseline_uncertainty=self.baseline_uncertainty.fetch(),
            baseline_name=self.baseline_name.fetch(),
            active_thresh=self.active_thresh.fetch(),
            active_out=self.active_output.fetch(),
            parse_json=self.parse_json.fetch(),
        )


class InputFFRotations(InputForceField):
    default_help = """
Wraps around another forcefield to evaluate it over one or more 
rotated copies of the physical system. This is useful when 
interacting with models that are not exactly invariant/covariant 
with respect to rigid rotations.
Besides the parameters defining how averaging is to be performed
(using an integration grid, and/or randomizing the orientation at
each step) the <ffrotations> should contain either a <ffsocket>
or a <ffdirect> block that computes the "base" model. Note that
this forcefield should be given a separate `name`, but that you
cannot access this "inner" forcefield from other parts of the 
input file.
"""
    default_label = "FFROTATIONS"

    fields = copy(InputForceField.fields)

    fields["random"] = (
        InputValue,
        {
            "dtype": bool,
            "default": False,
            "help": """Applies a random rotation at each evaluation. """,
        },
    )

    fields["grid_order"] = (
        InputValue,
        {
            "dtype": int,
            "default": 1,
            "help": """Sums over a grid of rotations of the given order.
            Note that the number of rotations increases rapidly with the order:
            e.g. for a legendre grid
            '1' leads to a single rotation, '2' to 18, '3' to 75 rotations, while 
            for a lebedev grid '3' contains 18 rotations, '5' 70 rotations and so on.
            """,
        },
    )

    fields["grid_mode"] = (
        InputValue,
        {
            "dtype": str,
            "options": ["legendre", "lebedev"],
            "default": "legendre",
            "help": """Defines the type of integration grid. 
            Lebedev grids are usually more efficient in integrating. 
            """,
        },
    )

    fields["inversion"] = (
        InputValue,
        {
            "dtype": bool,
            "default": False,
            "help": """Always applies the improper version of each rotation in the
            grid (or the randomly-sampled rotation). Doubles the evaluations and makes
            the model exactly equivariant to inversion.""",
        },
    )

    fields["prng"] = (
        InputRandom,
        {
            "help": InputRandom.default_help,
            "default": input_default(factory=Random),
        },
    )

    fields["ffsocket"] = (
        InputFFSocket,
        {
            "help": InputFFSocket.default_help,
            "default": input_default(factory=FFSocket, kwargs={"name": "__DUMMY__"}),
        },
    )

    fields["ffdirect"] = (
        InputFFDirect,
        {
            "help": InputFFDirect.default_help,
            "default": input_default(factory=FFDirect, kwargs={"name": "__DUMMY__"}),
        },
    )

    def store(self, ff):
        """Store the base and rotation quadrature parameters"""

        super(InputFFRotations, self).store(ff)
        self.inversion.store(ff.inversion)
        self.grid_order.store(ff.grid_order)
        self.grid_mode.store(ff.grid_mode)
        self.random.store(ff.random)
        self.prng.store(ff.prng)
        self.ffsocket.store(ff.ffsocket)
        self.ffdirect.store(ff.ffdirect)

    def fetch(self):
        """Fetches all of the FF objects"""

        return FFRotations(
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            active=self.activelist.fetch(),
            threaded=self.threaded.fetch(),
            prng=self.prng.fetch(),
            ffsocket=self.ffsocket.fetch(),
            ffdirect=self.ffdirect.fetch(),
            grid_order=self.grid_order.fetch(),
            grid_mode=self.grid_mode.fetch(),
            random=self.random.fetch(),
            inversion=self.inversion.fetch(),
        )


class InputFFCavPhSocket(InputFFSocket):
    default_help = """A cavity molecular dynamics driver for vibraitonal strong coupling. 
                      In the current implementation, only a single cavity mode polarized along the x and y directions is coupled to the molecules.
                      Check https://doi.org/10.1073/pnas.2009272117 and also examples/lammps/h2o-cavmd/ for details.
                   """
    default_label = "FFCAVPHSOCKET"

    fields = {
        "charge_array": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "The partial charges of all the atoms, in the format [Q1, Q2, ... ].",
                "dimension": "length",
            },
        ),
        "apply_photon": (
            InputValue,
            {
                "dtype": bool,
                "default": False,
                "help": "Determines if additional photonic degrees of freedom is included or not.",
            },
        ),
        "E0": (
            InputValue,
            {
                "dtype": float,
                "default": 0.0,
                "help": "The value of varepsilon (effective light-matter coupling strength) in atomic units.",
            },
        ),
        "omega_c": (
            InputValue,
            {
                "dtype": float,
                "default": 0.01,
                "help": "This gives the cavity photon frequency at normal incidence.",
                "dimension": "frequency",
            },
        ),
        "ph_rep": (
            InputValue,
            {
                "dtype": str,
                "default": "loose",
                "options": [
                    "loose",
                    "dense",
                ],
                "help": """In the current implementation, two energy-degenerate photon modes polarized along x and y directions
                are coupled to the molecular system. If 'loose', the cavity photons polarized along the x, y directions are represented by two 'L' atoms; 
                the x dimension of the first 'L' atom is coupled to the molecules, and the y dimension of the second 'L' atom is coupled to the molecules.
                If 'dense', the cavity photons polarized along the x, y directions are represented by one 'L' atom; 
                the x and y dimensions of this 'L' atom are coupled to the molecules.""",
            },
        ),
    }

    fields.update(InputFFSocket.fields)

    attribs = {}
    attribs.update(InputFFSocket.attribs)

    def store(self, ff):
        """Takes a ForceField instance and stores a minimal representation of it.

        Args:
           ff: A ForceField object with a FFCavPhSocket forcemodel object.
        """

        if not type(ff) is FFCavPhSocket:
            raise TypeError(
                "The type " + type(ff).__name__ + " is not a valid socket forcefield"
            )
        super(InputFFCavPhSocket, self).store(ff)

        self.charge_array.store(ff.charge_array)
        self.apply_photon.store(ff.apply_photon)
        self.E0.store(ff.E0)
        self.omega_c.store(ff.omega_c)
        self.ph_rep.store(ff.ph_rep)

    def fetch(self):
        """Creates a ForceSocket object.

        Returns:
           A ForceSocket object with the correct socket parameters.
        """

        if self.threaded.fetch() == False:
            raise ValueError("FFCavPhFPSockets cannot poll without threaded mode.")

        # just use threaded throughout
        return FFCavPhSocket(
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            offset=self.offset.fetch(),
            dopbc=self.pbc.fetch(),
            active=self.activelist.fetch(),
            threaded=self.threaded.fetch(),
            interface=InterfaceSocket(
                address=self.address.fetch(),
                port=self.port.fetch(),
                slots=self.slots.fetch(),
                mode=self.mode.fetch(),
                timeout=self.timeout.fetch(),
                match_mode=self.matching.fetch(),
                exit_on_disconnect=self.exit_on_disconnect.fetch(),
            ),
            charge_array=self.charge_array.fetch(),
            apply_photon=self.apply_photon.fetch(),
            E0=self.E0.fetch(),
            omega_c=self.omega_c.fetch(),
            ph_rep=self.ph_rep.fetch(),
        )
