"""Replay of the MD steps as CUDA graphs (see ipi.utils.cudagraphs).

The step of the integrator of a system is recorded as CUDA graphs, and
replayed at the following steps. Force evaluations are not recorded, as they
involve communication with forcefields and drivers that may synchronize with
the device: the recording is split in segments at each force evaluation,
which is done eagerly between the replays of the segments. For direct
forcefields whose driver works on the arrays of the backend, only the calls
to the driver are eager: the preparation of the requests and the collection
of the results are recorded, and the results of the driver are copied into
those it returned while the step was recorded.

Only the dynamics that are known to work in this way are replayed (see
`unsupported`), and a `CUDAGraphError` stops the simulation if a host value
read by a graph changes.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

from functools import partial

from ipi.utils.cudagraphs import CUDAGraphError, GraphReplay
from ipi.utils.messages import verbosity, warning

__all__ = ["replay_steps"]


def unsupported(simulation, system):
    """Returns why the steps of `system` cannot be replayed as CUDA graphs, or
    None if they can."""

    from ipi.engine.motion.dynamics import Dynamics, NVEIntegrator, NVTIntegrator
    from ipi.engine.smotion import Smotion
    from ipi.engine.thermostats import ThermoLangevin, ThermoPILE_L, Thermostat

    motion = system.motion
    if type(motion) is not Dynamics:
        return f"the motion class '{type(motion).__name__}'"
    if type(motion.integrator) not in (NVEIntegrator, NVTIntegrator):
        return f"the ensemble '{motion.enstype}'"
    if type(motion.thermostat) not in (Thermostat, ThermoLangevin, ThermoPILE_L):
        return f"the thermostat '{type(motion.thermostat).__name__}'"
    if len(motion.fixatoms_dof) > 0:
        return "fixed atoms"
    if len(system.nm.open_paths) > 0 or len(system.nm.bosons) > 0:
        return "open paths or bosons"
    if simulation.smotion is not None and type(simulation.smotion) is not Smotion:
        return "a motion of multiple systems"
    if simulation.threading and len(simulation.syslist) > 1:
        return "stepping systems in threads"
    return None


def replay_steps(simulation):
    """Replays the steps of the systems of `simulation` as CUDA graphs, for the
    systems that support it."""

    for system in simulation.syslist:
        reason = unsupported(simulation, system)
        if reason is None:
            GraphStep(system)
        else:
            warning(
                f"The steps are not replayed as CUDA graphs: {reason} is not supported",
                verbosity.low,
            )


def _eager_driver(ff):
    """Whether the requests to `ff` are recorded, with only the calls to its
    driver done eagerly."""

    from ipi.engine.forcefields import FFDirect

    return type(ff) is FFDirect and ff.backend_arrays and not ff.threaded


class GraphStep(GraphReplay):
    """Replaces the step of the integrator of a system: runs the first steps
    eagerly, records the next one as CUDA graphs, and replays them for all the
    following steps."""

    def __init__(self, system):
        integrator = system.motion.integrator
        super().__init__(integrator.step, system.prng.torch_generator)
        integrator.step = self
        self.prng = system.prng

        # the force components, which may be computed eagerly, of each bead
        self.components = {}
        for forces in (system.forces, system.ensemble.bias):
            for component in getattr(forces, "mforces", []):
                # scaled components are computed from a base component
                component = getattr(component, "bf", component)
                for bead in component._forces:
                    self.components[bead] = component

    def __call__(self, step=None):
        if (
            self.segments is not None
            and self.prng.torch_generator is not self.generator
        ):
            raise CUDAGraphError("The random number generator has been replaced")
        super().__call__(step)

    def split(self, bead):
        """Interrupts the recording to compute the forces on `bead` eagerly,
        unless they are already computed, or only the calls to the driver of
        its forcefield are eager."""

        if bead._ufvx.tainted() and not _eager_driver(bead.ff):
            self.interrupt(partial(_compute, self.components[bead]))


def _compute(component):
    """Computes the forces, potentials and virials of a force component."""

    component._f.get()
    component._pots.get()
    component._virs.get()
