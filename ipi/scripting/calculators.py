"""Python-side calculators used by the scripting API.

Currently provides a minimal ASE-compatible calculator used internally
by :class:`ipi.scripting.interactive.InteractiveSimulation` when
returning simulation snapshots as ASE ``Atoms`` objects. Future
additions (e.g. a direct metatomic-based calculator that does not go
through the ASE compatibility layer) will live here as well.
"""

from ipi.utils.io.backends.io_ase import _asecheck


DummyASECalculator = None


def _define_calculator():
    global DummyASECalculator
    _asecheck()
    from ase.calculators.calculator import Calculator, all_changes

    if DummyASECalculator is None:

        class DummyASECalculator(Calculator):
            implemented_properties = ["energy", "forces"]

            def __init__(self, energy=None, forces=None, **kwargs):
                super().__init__(**kwargs)
                self._energy = energy
                self._forces = forces

            def calculate(
                self,
                atoms=None,
                properties=["energy", "forces"],
                system_changes=all_changes,
            ):
                super().calculate(atoms, properties, system_changes)
                self.results = {}
                if "energy" in properties:
                    self.results["energy"] = self._energy
                if "forces" in properties:
                    self.results["forces"] = self._forces
