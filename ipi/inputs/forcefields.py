"""Creates objects that deal with the evaluation of interactions."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from copy import copy
import numpy as np

from ipi.engine.forcefields import (
    ForceField,
    FFSocket,
    FFLennardJones,
    FFDebye,
    FFPlumed,
    FFYaff,
    FFsGDML,
    FFCommittee,
)
from ipi.interfaces.sockets import InterfaceSocket
import ipi.engine.initializer
from ipi.inputs.initializer import *
from ipi.utils.inputvalue import *

__all__ = [
    "InputFFSocket",
    "InputFFLennardJones",
    "InputFFDebye",
    "InputFFPlumed",
    "InputFFYaff",
    "InputFFsGDML",
    "InputFFCommittee",
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
          code.

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
                "default": True,
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
                "default": 0.01,
                "help": "The number of seconds the polling thread will wait between exhamining the list of requests.",
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
                "help": "List with indexes of the atoms that this socket is taking care of.    Default: all (corresponding to -1)",
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
        self.parameters.store(ff.pars)
        self.pbc.store(ff.dopbc)
        self.activelist.store(ff.active)
        self.threaded.store(ff.threaded)

    def fetch(self):
        """Creates a ForceField object.

        Returns:
           A ForceField object.
        """

        super(InputForceField, self).fetch()

        return ForceField(
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
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
                "options": ["auto", "any"],
                "default": "auto",
                "help": "Specifies whether requests should be dispatched to any client, or automatically matched to the same client when possible [auto].",
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

        if not type(ff) is FFSocket:
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
        return FFSocket(
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
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


class InputFFLennardJones(InputForceField):

    attribs = {}
    attribs.update(InputForceField.attribs)

    default_help = """Simple, internal LJ evaluator without cutoff, neighbour lists or minimal image convention.
                   Expects standard LJ parameters, e.g. { eps: 0.1, sigma: 1.0 }. """
    default_label = "FFLJ"

    def store(self, ff):
        super(InputFFLennardJones, self).store(ff)

    def fetch(self):
        super(InputFFLennardJones, self).fetch()

        return FFLennardJones(
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            dopbc=self.pbc.fetch(),
            threaded=self.threaded.fetch(),
        )

        if self.slots.fetch() < 1 or self.slots.fetch() > 5:
            raise ValueError(
                "Slot number " + str(self.slots.fetch()) + " out of acceptable range."
            )
        if self.latency.fetch() < 0:
            raise ValueError("Negative latency parameter specified.")
        if self.timeout.fetch() < 0.0:
            raise ValueError("Negative timeout parameter specified.")


class InputFFDebye(InputForceField):

    fields = {
        "hessian": (
            InputArray,
            {
                "dtype": float,
                "default": input_default(factory=np.zeros, args=(0,)),
                "help": "Specifies the Hessian of the harmonic potential. "
                "Default units are atomic. Units can be specified only by xml attribute. "
                "Implemented options are: 'atomic_unit', 'ev/ang\^2'",
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
            dopbc=self.pbc.fetch(),
            threaded=self.threaded.fetch(),
        )


class InputFFPlumed(InputForceField):

    fields = {
        "init_file": (
            InputInitFile,
            {
                "default": input_default(
                    factory=ipi.engine.initializer.InitFile, kwargs={"mode": "xyz"}
                ),
                "help": "This describes the location to read the reference structure file from.",
            },
        ),
        "plumeddat": (
            InputValue,
            {"dtype": str, "default": "plumed.dat", "help": "The PLUMED input file"},
        ),
        "plumedstep": (
            InputValue,
            {
                "dtype": int,
                "default": 0,
                "help": "The current step counter for PLUMED calls",
            },
        ),
    }

    attribs = {}

    attribs.update(InputForceField.attribs)
    fields.update(InputForceField.fields)

    default_help = """ Direct PLUMED interface """
    default_label = "FFPLUMED"

    def store(self, ff):
        super(InputFFPlumed, self).store(ff)
        self.plumeddat.store(ff.plumeddat)
        # pstep = ff.plumedstep
        # if pstep > 0: pstep -= 1 # roll back plumed step before writing a restart
        # self.plumedstep.store(pstep)
        self.plumedstep.store(ff.plumedstep)
        self.init_file.store(ff.init_file)

    def fetch(self):
        super(InputFFPlumed, self).fetch()

        return FFPlumed(
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            dopbc=self.pbc.fetch(),
            threaded=self.threaded.fetch(),
            plumeddat=self.plumeddat.fetch(),
            plumedstep=self.plumedstep.fetch(),
            init_file=self.init_file.fetch(),
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
                      Also contains options to use it for uncertainty estimation and for
                      active learning in a ML context, based on a committee model.
                      Implements the approaches discussed in DOI: 10.1063/5.0036522.
                      """
    default_label = "FFCOMMITTEE"

    dynamic = {
        "ffsocket": (InputFFSocket, {"help": InputFFSocket.default_help}),
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

    def store(self, ff):
        """ Store all the sub-forcefields """

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

        for (_ii, _obj) in enumerate(_fflist):
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
        """ Fetches all of the FF objects """

        super(InputFFCommittee, self).fetch()

        fflist = []
        for (k, v) in self.extra:
            fflist.append(v.fetch())

        # TODO: will actually need to create a FF object here!
        return FFCommittee(
            pars=self.parameters.fetch(),
            name=self.name.fetch(),
            latency=self.latency.fetch(),
            dopbc=self.pbc.fetch(),
            fflist=fflist,
            ffweights=self.weights.fetch(),
            alpha=self.alpha.fetch(),
            baseline_uncertainty=self.baseline_uncertainty.fetch(),
            baseline_name=self.baseline_name.fetch(),
            active_thresh=self.active_thresh.fetch(),
            active_out=self.active_output.fetch(),
        )
