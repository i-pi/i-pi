"""Inputs and analytic drivers for small reference simulations."""

import json
import xml.etree.ElementTree as ET

import numpy as np


def harmonic_input(nbeads=1, thermostat=None, threading=False, timestep=0.02):
    """Two atoms (masses 2 and 3) in U = 0.25 |q|^2, in atomic units."""
    q = np.tile([0.2, -0.3, 0.1, -0.4, 0.2, 0.3], (nbeads, 1))
    q[:, 0] += np.arange(nbeads) * 0.01
    p = np.tile([0.1, 0.2, -0.1, 0.3, -0.2, 0.1], (nbeads, 1))
    bath = (
        ""
        if thermostat is None
        else f'<thermostat mode="{thermostat}"><tau>0.5</tau></thermostat>'
    )
    mode = "nve" if thermostat is None else "nvt"
    return f"""<simulation verbosity="quiet" threading="{threading}" floatformat="%24.16e">
      <total_steps>100</total_steps><prng><seed>317</seed></prng>
      <ffdirect name="harmonic"><pes>harmonic</pes><parameters>{{k1: 0.5}}</parameters></ffdirect>
      <system><forces><force forcefield="harmonic"/></forces>
        <ensemble><temperature units="kelvin">300</temperature></ensemble>
        <motion mode="dynamics"><fixcom>False</fixcom>
          <dynamics mode="{mode}"><timestep>{timestep}</timestep>{bath}</dynamics>
        </motion>
        <beads natoms="2" nbeads="{nbeads}">
          <q shape="({nbeads},6)">{q.ravel().tolist()}</q>
          <p shape="({nbeads},6)">{p.ravel().tolist()}</p>
          <m shape="(2)">[2, 3]</m><names shape="(2)">[H, O]</names>
        </beads><cell shape="(3,3)">[20,0,0,0,20,0,0,0,20]</cell>
      </system></simulation>"""


def motion_input(mode, options, nbeads=1):
    """Harmonic input tree whose motion is replaced by `mode` with `options`."""
    root = ET.fromstring(harmonic_input(nbeads=nbeads))
    motion = root.find("system/motion")
    motion.clear()
    motion.set("mode", mode)
    ET.SubElement(motion, "fixcom").text = "False"
    motion.append(options)
    return root


def pressure_input(barostat, nbeads=1, constant_volume=False):
    """Harmonic input with a Langevin-thermostatted barostat of the given mode."""
    root = ET.fromstring(
        harmonic_input(nbeads=nbeads, thermostat="langevin", timestep=0.001)
    )
    ensemble = root.find("system/ensemble")
    ET.SubElement(ensemble, "pressure").text = ".000001"
    dynamics = root.find("system/motion/dynamics")
    dynamics.set("mode", "nst" if barostat == "anisotropic" else "npt")
    if barostat == "anisotropic":
        ET.SubElement(ensemble, "stress", shape="(3,3)").text = str(
            (np.eye(3) * 0.000001).ravel().tolist()
        )
    baro = ET.SubElement(dynamics, "barostat", mode=barostat)
    ET.SubElement(baro, "tau").text = "100"
    thermo = ET.SubElement(baro, "thermostat", mode="langevin")
    ET.SubElement(thermo, "tau").text = "100"
    if barostat == "anisotropic":
        ET.SubElement(baro, "h0", shape="(3,3)").text = str(
            (np.eye(3) * 20).ravel().tolist()
        )
    if constant_volume:
        ET.SubElement(baro, "vol_constraint").text = "True"
    return ET.tostring(root, encoding="unicode")


class PrescribedPolynomial:
    """Driver for U = sum(a x^2 + b x^4) over all Cartesian coordinates."""

    def __init__(self, quadratic, quartic):
        self.a = np.array(quadratic)
        self.b = np.array(quartic)

    def __call__(self, cell, positions):
        x = positions.reshape(-1)
        energy = np.sum(self.a * x**2 + self.b * x**4)
        force = -2 * self.a * x - 4 * self.b * x**3
        virial = positions.T @ force.reshape(-1, 3)
        extra = json.dumps({"energy": float(energy)})
        return energy, force.reshape(positions.shape), virial, extra


# Cases of the tests of the backends and of the compiled steps (--jit), which
# compare trajectories: thermostats, barostats and dynamics on a box of hydrogen
# atoms, written by `write_input` (or by the functions of INPUTS).

LANGEVIN = "<thermostat mode='langevin'><tau units='femtosecond'>100</tau></thermostat>"
SVR = "<thermostat mode='svr'><tau units='femtosecond'>100</tau></thermostat>"
PILE_L = "<thermostat mode='pile_l'><tau units='femtosecond'>100</tau></thermostat>"
PILE_G = "<thermostat mode='pile_g'><tau units='femtosecond'>100</tau></thermostat>"
GLE = """<thermostat mode='gle'><A shape='(3,3)'>
  [ 1e-3, 2e-3, 1e-4, -2e-3, 3e-3, 1e-4, -1e-4, -1e-4, 2e-3 ]</A></thermostat>"""
FFL = "<thermostat mode='ffl'><tau units='femtosecond'>100</tau></thermostat>"
CL = """<thermostat mode='cl'><tau units='femtosecond'>100</tau>
  <intau units='femtosecond'>200</intau></thermostat>"""
MULTI = f"<thermostat mode='multi'>{LANGEVIN}{SVR}</thermostat>"


def barostat_xml(mode, thermostat=LANGEVIN, extra=""):
    return f"""<barostat mode='{mode}'><tau units='femtosecond'>200</tau>
  {thermostat}{extra}</barostat>"""


ISOTROPIC = barostat_xml("isotropic")
FLEXIBLE = barostat_xml("flexible")
# a second forcefield, for multiple time steps, contraction and bias
SECOND = {
    "ff": "<ffdirect name='ff2'><pes>harmonic</pes>"
    "<parameters>{k1: 5e-4}</parameters></ffdirect>"
}
MTS = {
    **SECOND,
    "forces": "<force forcefield='ff'><mts_weights>[0, 1]</mts_weights></force>"
    "<force forcefield='ff2'><mts_weights>[1, 0]</mts_weights></force>",
    "dynamics": "<nmts>[1, 2]</nmts>",
}

# (number of beads, ensemble, thermostat, barostat[, extra input sections])
CASES = {
    "nve": (1, "nve", "", ""),
    "nvt_langevin": (1, "nvt", LANGEVIN, ""),
    "nvt_svr": (1, "nvt", SVR, ""),
    "nvt_gle": (1, "nvt", GLE, ""),
    "nvt_ffl": (1, "nvt", FFL, ""),
    "nvt_cl": (1, "nvt", CL, ""),
    "nvt_multi": (1, "nvt", MULTI, ""),
    "nvt_mts_bias": (
        1,
        "nvt",
        LANGEVIN,
        "",
        {**MTS, "ensemble": "<bias><force forcefield='ff2'/></bias>"},
    ),
    "npt_isotropic": (1, "npt", LANGEVIN, ISOTROPIC),
    "npt_flexible": (1, "npt", LANGEVIN, FLEXIBLE),
    "npt_svr_barostat": (1, "npt", LANGEVIN, barostat_xml("flexible", SVR)),
    "npt_gle_barostat": (1, "npt", LANGEVIN, barostat_xml("flexible", GLE)),
    "nst_anisotropic": (
        1,
        "nst",
        LANGEVIN,
        barostat_xml(
            "anisotropic", extra="<h0 shape='(3,3)'>[10, 0, 0, 0, 10, 0, 0, 0, 10]</h0>"
        ),
        {
            "ensemble": "<stress shape='(3,3)'>[1e-6, 0, 0, 0, 2e-6, 0, 0, 0, 1e-6]</stress>"
        },
    ),
    "pimd_nvt": (4, "nvt", PILE_L, ""),
    "pimd_pile_g": (4, "nvt", PILE_G, ""),
    "pimd_npt": (4, "npt", PILE_L, ISOTROPIC),
    "pimd_nvt_cc": (4, "nvt-cc", PILE_L, ""),
    "pimd_fixatoms": (
        4,
        "nvt",
        PILE_G,
        "",
        {"motion": "<fixcom>False</fixcom><fixatoms>[0, 3]</fixatoms>"},
    ),
    "pimd_sc": (4, "sc", PILE_L, ""),
    "pimd_sc_forward": (
        4,
        "sc",
        PILE_L,
        "",
        {"forces": "<force forcefield='ff' fd_epsilon='1e-3'/>"},
    ),
    "pimd_scnpt": (4, "scnpt", PILE_L, barostat_xml("sc-isotropic")),
    "pimd_contracted": (
        4,
        "nvt",
        PILE_L,
        "",
        {
            **SECOND,
            "forces": "<force forcefield='ff' nbeads='2'/><force forcefield='ff2'/>",
        },
    ),
}

HARMONIC = (
    "<ffdirect name='ff'><pes>harmonic</pes>"
    "<parameters>{k1: 1e-3}</parameters></ffdirect>",
    "<force forcefield='ff'/>",
)
# a force component with zero weight is never evaluated
NO_FORCES = (
    "<fflj name='ff' pbc='False'><parameters>{eps: 0.0, sigma: 1.0}</parameters></fflj>",
    "<force forcefield='ff' weight='0.0'/>",
)


def write_input(case, forcefield, natoms=8):
    """Writes `init.xyz` in the current folder and returns the XML input."""

    nbeads, mode, thermostat, barostat, *extra = CASES[case]
    extra = extra[0] if extra else {}
    box = (natoms / 0.05) ** (1.0 / 3.0)
    pos = np.random.default_rng(12345).random((natoms, 3)) * box
    with open("init.xyz", "w") as f:
        f.write(f"{natoms}\n# CELL(abcABC): {box} {box} {box} 90 90 90 ")
        f.write("cell{angstrom} positions{angstrom}\n")
        for x, y, z in pos:
            f.write(f"H {x:.8f} {y:.8f} {z:.8f}\n")
    pressure = "<pressure units='bar'>1</pressure>" if barostat else ""
    ffxml, forcexml = forcefield
    forcexml = extra.get("forces", forcexml)
    return f"""
<simulation verbosity='quiet' threading='false' safe_stride='100000000'>
  {ffxml}{extra.get("ff", "")}
  <total_steps>100000000</total_steps>
  <output prefix='simulation'/>
  <prng><seed>12345</seed></prng>
  <system>
    <initialize nbeads='{nbeads}'>
      <file mode='xyz'>init.xyz</file>
      <velocities mode='thermal' units='kelvin'>300</velocities>
    </initialize>
    <forces>{forcexml}</forces>
    <ensemble>
      <temperature units='kelvin'>300</temperature>{pressure}{extra.get("ensemble", "")}
    </ensemble>
    <motion mode='dynamics'>{extra.get("motion", "")}
      <dynamics mode='{mode}'>
        <timestep units='femtosecond'>0.5</timestep>{thermostat}{barostat}
        {extra.get("dynamics", "")}
      </dynamics>
    </motion>
  </system>
</simulation>
"""


def mts_input(forcefield):
    """PIMD with a forcefield and the dummy potential applied at two levels of
    multiple time stepping."""

    xml = write_input("pimd_nvt", forcefield)
    return (
        xml.replace(
            forcefield[0],
            forcefield[0] + "<ffdirect name='ff2'><pes>dummy</pes></ffdirect>",
        )
        .replace(
            "<force forcefield='ff'/>",
            "<force forcefield='ff'><mts_weights>[1,0]</mts_weights></force>"
            "<force forcefield='ff2'><mts_weights>[0,1]</mts_weights></force>",
        )
        .replace(
            "<timestep units='femtosecond'>0.5</timestep>",
            "<timestep units='femtosecond'>0.5</timestep><nmts>[1,3]</nmts>",
        )
    )


def nst_input(forcefield):
    """Classical dynamics at constant stress, with the anisotropic barostat."""

    box = (8 / 0.05) ** (1.0 / 3.0)  # the box of write_input, in angstrom
    return (
        write_input("npt_isotropic", forcefield)
        .replace(
            "<pressure units='bar'>1</pressure>",
            "<stress units='bar'>[1, 0, 0, 0, 1, 0, 0, 0, 1]</stress>",
        )
        .replace("mode='npt'", "mode='nst'")
        .replace(
            "<barostat mode='isotropic'>",
            "<barostat mode='anisotropic'><h0 units='angstrom'>"
            f"[{box}, 0, 0, 0, {box}, 0, 0, 0, {box}]</h0>",
        )
    )


INPUTS = {"mts": mts_input, "nst": nst_input}

DUMMY = ("<ffdirect name='ff'><pes>dummy</pes></ffdirect>", "<force forcefield='ff'/>")
DUMMY_BATCHED = (
    "<ffdirect name='ff'><pes>dummy</pes><batch_size>4</batch_size></ffdirect>",
    "<force forcefield='ff'/>",
)
# the properties compared between compiled and eager steps
PROPERTIES = ["kinetic_md", "potential", "conserved", "temperature"]
# compiled steps match eager ones to rounding errors, as the compiled functions
# may sum in a different order
TOLERANCE = dict(rtol=1e-12, atol=1e-14)


def case_input(case, forcefield):
    """The input of a case, from CASES or INPUTS, with `forcefield`."""

    return INPUTS[case](forcefield) if case in INPUTS else write_input(case, forcefield)
