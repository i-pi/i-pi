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
