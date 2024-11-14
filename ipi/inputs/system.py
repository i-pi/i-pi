"""Deals with creating a representation of a system."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import ipi.engine.system
from ipi.utils.depend import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *
from ipi.utils.prng import *
from ipi.utils.io import *
from ipi.utils.io.inputs.io_xml import *
from ipi.inputs.forces import InputForces
from ipi.inputs.beads import InputBeads
from ipi.inputs.cell import InputCell
from ipi.inputs.ensembles import InputEnsemble
from ipi.inputs.motion import InputMotion
from ipi.inputs.normalmodes import InputNormalModes
from ipi.engine.normalmodes import NormalModes
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.engine.ensembles import Ensemble
from ipi.engine.motion import Motion
from ipi.inputs.initializer import InputInitializer
from ipi.engine.initializer import Initializer

__all__ = ["InputSystem", "InputSysTemplate"]


class InputSysTemplate(Input):
    """A template class to generate multiple systems with varying parameters.

    Generates a series of systems by automatically filling up labels in a template input structure

    Fields:
        template: The text that corresponds to the system field template
        labels: List of strings identifying substitution fields in the template
        instance: List of strings that should be used to substitute the labels
    """

    fields = {
        "template": (
            InputRaw,
            {
                "help": """ A string that will be read verbatim containing the model for a system to be generated""",
                "dtype": str,
            },
        ),
        "labels": (
            InputArray,
            {
                "help": """ A list of strings that should be substituted in the template to create multiple systems """,
                "dtype": str,
            },
        ),
    }
    dynamic = {
        "instance": (
            InputArray,
            {
                "help": """ A list of strings that should the labels creating one system instance """,
                "dtype": str,
            },
        )
    }

    def fetch(self):
        """Creates a series of physical system objects using a template and label substitutions.

        Returns:
            A list pf System objects of the appropriate type and with the appropriate
            properties, built by substituting placeholders in a template with the given values

        Raises:
            ValueError: Raised if the labels and the instance lists have mismatching lengths
        """

        super(InputSysTemplate, self).fetch()

        template = self.template.fetch()
        labels = self.labels.fetch()
        lsys = []
        for k, v in self.extra:
            if k == "instance":
                ins = v.fetch()
                sys = template
                if len(labels) != len(ins):
                    raise ValueError("Labels and instance length mismatch")
                for label, to_insert in sorted(zip(labels, ins), reverse=True):
                    # sort from longest to smallest label, to avoid replacing 'string' in 'string1'
                    sys = sys.replace(label, to_insert)
                print(" @inputsystemplate.fetch: Generating system from template:", sys)
                xsys = xml_parse_string(sys)  # parses the string to an XML object
                isys = InputSystem()
                isys.parse(
                    xsys.fields[0][1]
                )  # parses the XML object into an InputSystem
                lsys.append(
                    isys.fetch()
                )  # fetches the generated System and appends to the list

        return lsys


class InputSystem(Input):
    """Physical system input class.

    Handles generating the appropriate forcefield class from the xml input file,
    and generating the xml checkpoint tags and data from an instance of the
    object.

    Attributes:
       prefix: A string to prepend to the output file names for this system.

    Fields:
       forces: A restart force instance. Used as a model for all the replicas.
       ensemble: A restart ensemble instance.
       beads: A restart beads instance.
       normal_modes: Setup of normal mode integrator.
       cell: A restart cell instance.
       initialize: An array of strings giving all the quantities that should
          be output.
    """

    fields = {
        "initialize": (
            InputInitializer,
            {
                "help": InputInitializer.default_help,
                "default": input_default(factory=Initializer),
            },
        ),
        "forces": (InputForces, {"help": InputForces.default_help}),
        "ensemble": (
            InputEnsemble,
            {
                "help": InputEnsemble.default_help,
                "default": input_default(factory=Ensemble, kwargs={"temp": 1.0}),
            },
        ),
        "motion": (
            InputMotion,
            {
                "help": InputMotion.default_help,
                "default": input_default(factory=Motion),
            },
        ),
        "beads": (
            InputBeads,
            {
                "help": InputBeads.default_help,
                "default": input_default(
                    factory=Beads, kwargs={"natoms": 0, "nbeads": 0}
                ),
            },
        ),
        "normal_modes": (
            InputNormalModes,
            {
                "help": InputNormalModes.default_help,
                "default": input_default(factory=NormalModes, kwargs={"mode": "rpmd"}),
            },
        ),
        "cell": (
            InputCell,
            {"help": InputCell.default_help, "default": input_default(factory=Cell)},
        ),
    }

    attribs = {
        "prefix": (
            InputAttribute,
            {
                "help": "Prepend this string to output files generated for this system. ",
                "default": "",
                "dtype": str,
            },
        )
    }

    default_help = "This is the class which holds all the data which represents a single state of the system."
    default_label = "SYSTEM"

    def store(self, psys):
        """Takes a System instance and stores a minimal representation of it.

        Args:
           psys: A physical system object.
        """

        super(InputSystem, self).store()

        self.prefix.store(psys.prefix)
        self.forces.store(psys.fcomp)
        self.ensemble.store(psys.ensemble)
        self.motion.store(psys.motion)
        self.beads.store(psys.beads)
        self.normal_modes.store(psys.nm)
        self.cell.store(psys.cell)

    def fetch(self):
        """Creates a physical system object.

        Returns:
           A System object of the appropriate type and with the appropriate
           properties and other objects given the attributes of the
           InputSystem object.

        Raises:
           TypeError: Raised if one of the file types in the stride keyword
              is incorrect.
        """

        super(InputSystem, self).fetch()

        # this creates a simulation object which gathers all the little bits
        # TODO use named arguments since this list is a bit too long...
        rsys = ipi.engine.system.System(
            init=self.initialize.fetch(),
            beads=self.beads.fetch(),
            nm=self.normal_modes.fetch(),
            cell=self.cell.fetch(),
            fcomponents=self.forces.fetch(),
            ensemble=self.ensemble.fetch(),
            motion=self.motion.fetch(),
            prefix=self.prefix.fetch(),
        )

        return rsys
