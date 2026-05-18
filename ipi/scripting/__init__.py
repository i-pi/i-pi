"""User-facing Python API for i-PI.

Organization:

- :mod:`ipi.scripting.templates` — XML fragment builders for systems,
  forcefields and motion classes - to generate templates for simulations.
- :mod:`ipi.scripting.interactive` — :class:`InteractiveSimulation`, a
  Python-driven stepper over an i-PI simulation.
- :mod:`ipi.scripting.calculators` — Python-side calculators (currently
  an ASE-compatible dummy calculator used when returning snapshots as
  ASE ``Atoms``; home for other calculators).
- :mod:`ipi.scripting.parsing` — post-processing helpers: read i-PI
  output files into numpy arrays and trajectories into ASE ``Atoms``.

The most common symbols are re-exported here for convenience.
"""

from ipi.scripting.templates import (
    simulation_xml,
    motion_nvt_xml,
    forcefield_xml,
)
from ipi.scripting.interactive import InteractiveSimulation
from ipi.scripting.parsing import read_output, read_trajectory

__all__ = [
    "simulation_xml",
    "motion_nvt_xml",
    "forcefield_xml",
    "InteractiveSimulation",
    "read_output",
    "read_trajectory",
]
