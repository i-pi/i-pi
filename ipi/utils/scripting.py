import numpy as np
from ipi.utils.depend import dstrip
from ipi.utils.units import Elements, unit_to_internal, unit_to_user
from ipi.inputs.beads import InputBeads
from ipi.inputs.cell import InputCell
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.engine.simulation import Simulation
from ipi.utils.io.inputs.io_xml import (
    xml_node,
    xml_parse_string,
    xml_parse_file,
    xml_write,
    xml_edit,
)

from ipi.utils.io.backends.io_ase import ase, _asecheck

__all__ = ["setup_forcefield", "setup_nvt", "InteractiveSimulation"]


def setup_forcefield(
    name, mode="direct", parameters=None, pes=None, address=None, port=None
):
    if mode == "direct":
        xml_ff = f"""
<ffdirect name='{name}'>
<pes>{"dummy" if pes is None else pes}</pes>
{"" if parameters is None else f"<parameters>{parameters}</parameters>"}
</ffdirect>
"""
    elif mode == "unix" or mode == "inet":
        if address is None:
            raise ValueError("Must specify address for {mode} forcefield")
        if mode == "inet" and port is None:
            raise ValueError("Must specify port for {mode} forcefield")
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


def setup_nvt(structures, forcefield, temperature, timestep, thermostat=None):
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
    beads.q = q * unit_to_internal("length", "angstrom", 1.0)
    beads.names = np.array(structure.symbols)
    if "mass" in structure.arrays:
        beads.m = structure.arrays["mass"] * unit_to_internal("mass", "dalton", 1.0)
    else:
        beads.m = np.array([Elements.mass(n) for n in beads.names])

    input_beads = InputBeads()
    input_beads.store(beads)

    cell = Cell(
        h=np.array(structure.cell) * unit_to_internal("length", "angstrom", 1.0)
    )
    input_cell = InputCell()
    input_cell.store(cell)

    # gets ff name
    xml_ff = xml_parse_string(forcefield).fields[0][1]
    ff_name = xml_ff.attribs["name"]

    # sets up thermostat
    if thermostat is None:
        if nbeads == 1:
            # defaults to svr
            xml_thermostat = f"""
<thermostat mode='svr'>
    <tau units='femtosecond'> {10*timestep} </tau>
</thermostat>
"""
        else:
            # defaults to pile_g
            xml_thermostat = f"""
<thermostat mode='pile_g'>
    <tau units='femtosecond'> {10*timestep} </tau>
    <pile_lambda> 0.5 </pile_lambda>
</thermostat>
"""

    xml_string = f"""
<simulation>
{forcefield}
<system>
{input_beads.write("beads")}
{input_cell.write("cell")}
<initialize nbeads='{nbeads}'>
<velocities mode='thermal' units='kelvin'> {temperature} </velocities>
</initialize>
<ensemble>
<temperature units='kelvin'> {temperature} </temperature>
</ensemble>
<forces>
<force forcefield='{ff_name}'> </force>
</forces>
<motion mode='dynamics'>
<dynamics mode='nvt'>
<timestep units='femtosecond'> {timestep} </timestep>
{xml_thermostat}
</dynamics>
</motion>
</system>
</simulation>
    """

    return xml_string


class InteractiveSimulation:
    """A wrapper to `Simulation` that allows accessing
    properties in a "safe" way from Python."""

    def __init__(self, xml_input):
        self.simulation = Simulation.load_from_xml(xml_input)

    def run(self, steps=1):
        for istep in range(steps):
            self.simulation.run_step(self.simulation.step + istep)
        self.simulation.step += steps

    def properties(self, property):
        props = []
        for system in self.simulation.syslist:
            props.append(system.properties[property][0])
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
                    * unit_to_user("length", "angstrom", 1.0),
                    symbols=dstrip(system.beads.names),
                    cell=dstrip(system.cell.h)
                    * unit_to_user("length", "angstrom", 1.0),
                )
                struc.arrays["ipi_velocities"] = (
                    dstrip(system.beads.p[b]).reshape(-1, 3)
                    * unit_to_user("length", "angstrom", 1.0)
                    / unit_to_user("time", "femtosecond", 1.0)
                )
                struc.arrays["ipi_forces"] = (
                    dstrip(system.forces.f[b]).reshape(-1, 3)
                    * unit_to_user("energy", "electronvolt", 1.0)
                    / unit_to_user("length", "angstrom", 1.0)
                )
                struc.info["ipi_potential"] = dstrip(
                    system.forces.pots[b]
                ) * unit_to_user("energy", "electronvolt", 1.0)
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
                    "length", "angstrom", 1.0
                )
                if "ipi_velocities" in struc.arrays:
                    system.beads.p[b] = (
                        struc.arrays["ipi_velocities"].flatten()
                        / unit_to_user("length", "angstrom", 1.0)
                        * unit_to_user("time", "femtosecond", 1.0)
                    )
