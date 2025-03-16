import numpy as np
from ipi.utils.depend import dstrip
from ipi.utils.units import Elements, unit_to_internal, unit_to_user
from ipi.inputs.beads import InputBeads
from ipi.inputs.cell import InputCell
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.engine.simulation import Simulation
from ipi.utils.io.inputs.io_xml import xml_parse_string, xml_write, write_dict
from ipi.pes import __drivers__

from ipi.utils.io.backends.io_ase import ase, _asecheck

__all__ = [
    "simulation_xml",
    "motion_nvt_xml",
    "forcefield_xml",
    "InteractiveSimulation",
]


DEFAULT_OUTPUTS = """
  <output prefix='simulation'>
    <properties stride='2' filename='out'>  [ step, time{picosecond}, conserved{electronvolt}, temperature{kelvin}, potential{electronvolt} ] </properties>
    <trajectory filename='pos' stride='20' cell_units='angstrom'> positions{angstrom} </trajectory>
    <checkpoint stride='200'/>
  </output>
"""


def simulation_xml(
    structures,
    forcefield,
    motion,
    temperature=None,
    output=None,
    verbosity="quiet",
    safe_stride=20,
    prefix=None,
):
    """
    A helper function to generate an XML string for a basic i-PI
    simulation input.

    param structures: ase.Atoms|list(ase.Atoms) an Atoms object containing
        the initial structure for the system (or list of Atoms objects to
        initialize a path integral simulation with many beads). The
        structures should be in standardized form (cf. standard_cell from ASE).
    param forcefield: str An XML-formatted string describing the forcefield
        to be used in the simulation
    param motion: str An XML-formatted string describing the motion class
        to be used in the simulation
    temperature: Optional(float) A float specifying the temperature. If not
        specified, the temperature won't be set, and you should use a NVE
        motion class
    output: Optional(str) An  XML-formatted string describing the output
        A default list of outputs will be generated if not provided
    verbosity: Optional(str), default "quiet". The level of logging done
        to stdout
    safe_stride: Optional(int), default 20. How often the internal checkpoint
        data should be stored
    prefix: Optional(str) The prefix to be used for simulation outputs.
    """

    if type(structures) is list:
        nbeads = len(structures)
        natoms = len(structures[0])
        q = np.vstack([frame.positions.reshape(1, -1) for frame in structures])
        structure = structures[0]
    else:
        structure = structures
        nbeads = 1
        natoms = len(structure)
        q = structure.positions.flatten()

    beads = Beads(nbeads=nbeads, natoms=natoms)
    beads.q = q * unit_to_internal("length", "ase", 1.0)
    beads.names = np.array(structure.symbols)
    if "mass" in structure.arrays:
        beads.m = structure.arrays["mass"] * unit_to_internal("mass", "ase", 1.0)
    else:
        beads.m = np.array([Elements.mass(n) for n in beads.names])

    input_beads = InputBeads()
    input_beads.store(beads)

    cell = Cell(h=np.array(structure.cell).T * unit_to_internal("length", "ase", 1.0))
    input_cell = InputCell()
    input_cell.store(cell)

    # gets ff name
    xml_ff = xml_parse_string(forcefield).fields[0][1]
    ff_name = xml_ff.attribs["name"]

    # parses the outputs and overrides prefix
    if output is None:
        output = DEFAULT_OUTPUTS
    if prefix is not None:
        xml_output = xml_parse_string(output)
        if xml_output.fields[0][0] != "output":
            raise ValueError("the output parameter should be a valid 'output' block")
        xml_output.fields[0][1].attribs["prefix"] = prefix
        output = xml_write(xml_output)

    return f"""
<simulation verbosity='{verbosity}' safe_stride='{safe_stride}'>
{forcefield}
{output}
<system>
{input_beads.write("beads")}
{input_cell.write("cell")}
{(
    f"<initialize nbeads='{nbeads}'>"
    f"<velocities mode='thermal' units='ase'> {temperature} </velocities>"
    "</initialize>"
    f"<ensemble>"
    f"<temperature units='ase'> {temperature} </temperature>"
    "</ensemble>"
    if temperature is not None else ""
)}
<forces>
<force forcefield='{ff_name}'> </force>
</forces>
{motion}
</system>
</simulation>
"""


def forcefield_xml(
    name, mode="direct", parameters=None, pes=None, address=None, port=None
):
    """
    A helper function to generate an XML string for a forcefield block.
    """

    if mode == "direct":
        if pes is None:
            raise ValueError(f"Must specify 'pes' for {mode} forcefields")
        elif pes not in __drivers__.keys():
            raise ValueError(f"Invalid value {pes} for 'pes'")
        if parameters is None:
            parameters = ""
        elif type(parameters) is dict:
            parameters = f"<parameters>{write_dict(parameters)}</parameters>"
        else:
            parameters = f"<parameters>{parameters}</parameters>"

        xml_ff = f"""
<ffdirect name='{name}'>
<pes>{"dummy" if pes is None else pes}</pes>
{parameters}
</ffdirect>
"""
    elif mode == "unix" or mode == "inet":
        if address is None:
            raise ValueError("Must specify address for {mode} forcefields")
        if mode == "inet" and port is None:
            raise ValueError("Must specify port for {mode} forcefields")
        xml_ff = f"""
<ffsocket name='{name}'>
<address>{address}</address>
{f"<port>{port}</port>" if mode=="inet" else ""}
<latency> 1e-4 </latency>
</ffsocket>
"""
    else:
        raise ValueError(
            "Invalid forcefield mode, use ['direct', 'unix', 'inet'] or set up manually"
        )

    return xml_ff


def motion_nvt_xml(timestep, thermostat=None, path_integrals=False):
    """
    A helper function to generate an XML string for a MD simulation input.
    """

    # sets up thermostat
    if thermostat is None:
        if path_integrals:
            # defaults to pile_g
            xml_thermostat = f"""
<thermostat mode='pile_g'>
    <tau units='ase'> {10*timestep} </tau>
    <pile_lambda> 0.5 </pile_lambda>
</thermostat>
"""
        else:
            # defaults to svr
            xml_thermostat = f"""
<thermostat mode='svr'>
    <tau units='ase'> {10*timestep} </tau>
</thermostat>
"""

    return f"""
<motion mode="dynamics">
<dynamics mode="nvt">
<timestep units="ase"> {timestep} </timestep>
{xml_thermostat}
</dynamics>
</motion>
"""


class InteractiveSimulation:
    """A wrapper to `Simulation` that allows accessing
    properties and configurations in a "safe" way from Python.
    The object is initialized from an XML input, and supports
    a minimalistic API to step through the simulation, extract
    properties and structures, and set configurations interactively.

    :param xml_input: An open XML file, or an XML-formatted string
       containing the i-PI simulation that should be run.
    """

    def __init__(self, xml_input):
        self.simulation = Simulation.load_from_xml(xml_input)

    def run(self, steps=1, write_outputs=True):
        """Stepper through the simulation.
        Run the simulation interactively for the specified number of steps.

        :param steps: int, number of steps the simulation should be advanced by
        """

        print(self.simulation.step)
        self.simulation.tsteps = self.simulation.step + steps
        self.simulation.run(write_outputs=write_outputs)
        # saves the RESTART file
        self.simulation.chk.write(store=True)
        self.simulation.step += 1

    def run_step(self, step):
        self.simulation.run_step(step)

    def properties(self, property):
        """Fetches properties from the simulation state.

        :param property: the name of a property to be fetched from the simulation
        """
        props = []
        for system in self.simulation.syslist:
            value, dimension, unit = system.properties[property]
            props.append(value * unit_to_user(dimension, "ase", 1.0))
        if len(props) == 1:
            return props[0]
        else:
            return props

    def get_structures(self):
        _asecheck()

        sys_structures = []
        for system in self.simulation.syslist:
            structures = []
            for b in range(system.beads.nbeads):
                struc = ase.Atoms(
                    positions=dstrip(system.beads.q[b]).reshape(-1, 3)
                    * unit_to_user("length", "ase", 1.0),
                    symbols=dstrip(system.beads.names),
                    cell=dstrip(system.cell.h).T * unit_to_user("length", "ase", 1.0),
                )
                struc.arrays["ipi_velocities"] = dstrip(system.beads.p[b]).reshape(
                    -1, 3
                ) * unit_to_user("velocity", "ase", 1.0)
                struc.arrays["ipi_forces"] = dstrip(system.forces.f[b]).reshape(
                    -1, 3
                ) * unit_to_user("force", "ase", 1.0)
                struc.info["ipi_potential"] = dstrip(
                    system.forces.pots[b]
                ) * unit_to_user("energy", "ase", 1.0)
                structures.append(struc)
            if len(structures) == 1:  # flatten the structure if there's just one beads
                structures = structures[0]
            sys_structures.append(structures)
        if len(sys_structures) == 1:  # flatten the structure if there's just one system
            sys_structures = sys_structures[0]
        return sys_structures

    def set_structures(self, sys_structures):
        _asecheck()

        if len(self.simulation.syslist) == 1:
            sys_structures = [sys_structures]
        for system, structures in zip(self.simulation.syslist, sys_structures):
            if system.beads.nbeads == 1:
                structures = [structures]
            for b, struc in enumerate(structures):
                system.beads.q[b] = struc.positions.flatten() * unit_to_internal(
                    "length", "ase", 1.0
                )
                system.cell.h = struc.cell.T * unit_to_internal("length", "ase", 1.0)
                if "ipi_velocities" in struc.arrays:
                    system.beads.p[b] = struc.arrays[
                        "ipi_velocities"
                    ].flatten() * unit_to_internal("velocity", "ase", 1.0)
