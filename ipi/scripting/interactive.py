"""Interactive, Python-side driver for an i-PI simulation.

Exposes :class:`InteractiveSimulation`, a thin wrapper over
:class:`ipi.engine.simulation.Simulation` that lets user scripts step
through a simulation, inspect properties, swap structures, override
motion step functions, and inject custom properties.
"""

from copy import deepcopy

import numpy as np

from ipi.utils.depend import dstrip
from ipi.utils.units import unit_to_internal, unit_to_user, Constants
from ipi.engine.simulation import Simulation
from ipi.engine.motion import Dynamics
from ipi.utils.io.backends.io_ase import ase, _asecheck

from ipi.scripting import calculators

__all__ = ["InteractiveSimulation"]


class InteractiveSimulation(Simulation):
    """A wrapper to `Simulation` that allows accessing
    properties and configurations in a "safe" way from Python.
    The object is initialized from an XML input, and supports
    a minimalistic API to step through the simulation, extract
    properties and structures, and set configurations interactively.

    :param xml_input: An open XML file, or an XML-formatted string
       containing the i-PI simulation that should be run.
    : param custom_properties: Optional(dict) A dictionary of custom properties
       to be added to the simulation. The keys should be the property names,
       and the values should be either a callable (in which case the property
       is assumed to be of size 1), or a dictionary with the following keys:
         - 'func': the callable that computes the property, defined as a function
              taking `self` as the compulsory argument (where `self` is a reference
              to the `properties` object of a system, holding references to beads,
              forces, etc.)
         - 'size': Optional(int), the size of the property (default 1)
         - 'dimension': Optional(str), the dimension of the property
            (default 'undefined', could be 'length', 'energy', etc.)
         - 'help': Optional(str), a help string describing the property
    """

    def __init__(self, xml_input, custom_properties=None):
        # not-so-great way to initialize from the class method-generated
        # Simulation object. one should probably reconsider how Simulation
        # is initialized as this hack might lead to issues further down
        # the line, but for the moment all looks good.

        self._custom_properties = custom_properties or {}
        sim = Simulation.load_from_xml(
            xml_input, post_init_hook=self._add_custom_properties
        )
        self.__dict__.update(sim.__dict__)

    def _add_custom_properties(self, simulation):
        for name, property in self._custom_properties.items():
            # checks if properties have the right format
            if callable(property):
                property = {"func": property, "size": 1}
            elif not isinstance(property, dict):
                raise ValueError("Custom property should be a callable or a dict")
            if "func" not in property:
                raise ValueError("Custom property should contain a 'func' entry")
            if "size" not in property:
                # defaults to size 1
                property["size"] = 1
            if "dimension" not in property:
                property["dimension"] = "undefined"
            if "help" not in property:
                property["help"] = (
                    "Custom property, the devs didn't bother to add a description."
                )
            for s in simulation.syslist:
                # hooks the property function to the system instance
                sys_property = deepcopy(property)
                sys_property["func"] = sys_property["func"].__get__(
                    s.properties, s.properties.__class__
                )
                s.properties.property_dict[name] = sys_property

    def run(self, steps=1, write_outputs=True):
        """Stepper through the simulation.
        Run the simulation interactively for the specified number of steps.

        :param steps: int, number of steps the simulation should be advanced by
        """

        self.tsteps = self.step + steps
        super().run(write_outputs=write_outputs)
        # saves the RESTART file
        self.chk.write(store=True)
        self.step += 1

    def properties(self, property):
        """Fetches properties from the simulation state.

        :param property: the name of a property to be fetched from the simulation
        """
        props = []
        for system in self.syslist:
            value, dimension, unit = system.properties[property]
            props.append(value * unit_to_user(dimension, "ase", 1.0))
        if len(props) == 1:
            return props[0]
        else:
            return props

    def get_structures(self):
        _asecheck()
        calculators._define_calculator()

        sys_structures = []
        for system in self.syslist:
            structures = []
            for b in range(system.beads.nbeads):
                struc = ase.Atoms(
                    positions=dstrip(system.beads.q[b]).reshape(-1, 3)
                    * unit_to_user("length", "ase", 1.0),
                    symbols=dstrip(system.beads.names),
                    cell=dstrip(system.cell.h).T * unit_to_user("length", "ase", 1.0),
                )
                struc.set_momenta(
                    dstrip(system.beads.p[b]).reshape(-1, 3)
                    * unit_to_user("momentum", "ase", 1.0)
                )
                # set a dummy ASE calculator to be able to attach forces and
                # energies to the structure; these are updated at each step
                # when the user calls `get_structures()`, so they should always
                # reflect the current state of the simulation.
                struc.calc = calculators.DummyASECalculator(
                    energy=dstrip(system.forces.pots[b])
                    * unit_to_user("energy", "ase", 1.0),
                    forces=dstrip(system.forces.f[b]).reshape(-1, 3)
                    * unit_to_user("force", "ase", 1.0),
                )
                structures.append(struc)
            if len(structures) == 1:  # flatten the structure if there's just one beads
                structures = structures[0]
            sys_structures.append(structures)
        if len(sys_structures) == 1:  # flatten the structure if there's just one system
            sys_structures = sys_structures[0]
        return sys_structures

    def set_structures(self, sys_structures):
        _asecheck()

        if len(self.syslist) == 1:
            sys_structures = [sys_structures]
        for system, structures in zip(self.syslist, sys_structures):
            if system.beads.nbeads == 1:
                structures = [structures]
            for b, struc in enumerate(structures):
                system.beads.q[b] = struc.positions.flatten() * unit_to_internal(
                    "length", "ase", 1.0
                )
                system.cell.h = struc.cell.T * unit_to_internal("length", "ase", 1.0)
                if "momenta" in struc.arrays:
                    system.beads.p[b] = struc.arrays[
                        "momenta"
                    ].flatten() * unit_to_internal("momentum", "ase", 1.0)

    def thermalize_momenta(self, temperature_K=None):
        """Resets system momenta"""

        for s in self.syslist:
            if temperature_K is not None:
                temperature = unit_to_internal("energy", "kelvin", temperature_K)
            else:
                temperature = s.ensemble.temp

            # also correct for the change in kinetic energy
            s.ensemble.eens += s.nm.kin

            # initialize in the nm basis to handle PIMD mass scaling
            s.nm.pnm[:] = (
                self.prng.gvec((s.beads.nbeads, 3 * s.beads.natoms))
                * np.sqrt(dstrip(s.nm.dynm3))
                * np.sqrt(s.beads.nbeads * temperature * Constants.kb)
            )

            s.ensemble.eens -= s.nm.kin

            # apply momentum constraints
            if isinstance(s.motion, Dynamics):
                s.motion.integrator.pconstraints()

    def get_motion_step(self):
        """
        Fetches the step function(s) associated with the system motion classes.

        Returns a single function, or a list depending if there's just one system
        or many.
        """

        step_list = []
        for s in self.syslist:
            step_list.append(s.motion.step)

        if len(step_list) == 1:
            return step_list[0]
        else:
            return step_list

    def set_motion_step(self, custom_step):
        """
        Overrides the step function of the main motion class with a custom one.
        `custom_step` should be defined as

        ```
        def custom_step(self, step=None):
            do_something
        ```

        and `self` will be a reference to the motion class itself; depending on the
        original motion class, this might contain thermostats, barostats, etc.

        Call with a single function, or with a list if you have many systems and
        want to set a different function for each system (for whatever reason)
        """

        if not hasattr(custom_step, "__len__"):
            custom_step = [custom_step for s in self.syslist]

        for s, f in zip(self.syslist, custom_step):
            s.motion.__dict__["step"] = f.__get__(s.motion, s.motion.__class__)
