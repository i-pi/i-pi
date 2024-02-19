"""
The i-PI package.
"""

# expose some utility functions in a more direct way
from ipi.utils.parsing import read_output, read_trajectory

__all__ = [
    "clients",
    "engine",
    "inputs",
    "interfaces",
    "utils",
    "ipi_global_settings",
    "read_output",
    "read_trajectory",
]

ipi_global_settings = {"floatformat": "%16.8e"}
