"""
The i-PI module.
"""

# version is defined in ipi/_version.py, also read by the build (setup.cfg)
from ipi._version import __version__

ipi_global_settings = {"floatformat": "%16.8e"}

# expose some utility functions in a more direct way
from ipi.scripting import read_output, read_trajectory
from ipi.utils.setup import install_driver

__all__ = [
    "clients",
    "engine",
    "inputs",
    "interfaces",
    "utils",
    "pes",
    "ipi_global_settings",
    "install_driver",
    "read_output",
    "read_trajectory",
]

from ipi.engine.simulation import Simulation


class IPI:
    def __init__(self, xml_data):
        self._simulation = Simulation.load_from_xml(xml_data)
