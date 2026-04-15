"""Backward-compatibility shim: use :mod:`ipi.scripting` instead."""

import warnings

warnings.warn(
    "ipi.utils.scripting has moved to ipi.scripting; "
    "update your imports to `from ipi.scripting import ...`. "
    "This compatibility shim will be removed in a future release.",
    DeprecationWarning,
    stacklevel=2,
)

from ipi.scripting.templates import (  # noqa: F401
    simulation_xml,
    motion_nvt_xml,
    forcefield_xml,
)
from ipi.scripting.interactive import InteractiveSimulation  # noqa: F401

__all__ = [
    "simulation_xml",
    "motion_nvt_xml",
    "forcefield_xml",
    "InteractiveSimulation",
]
