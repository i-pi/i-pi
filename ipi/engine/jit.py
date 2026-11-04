"""Compilation of the MD steps (`--jit`): replay as CUDA graphs with torch (see
ipi.utils.cudagraphs), or as functions compiled by JAX (see ipi.utils.jaxjit).

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
`unsupported`), and a `JitError` stops the simulation if a host value read by
a graph changes.

With JAX, the whole step is compiled, including the requests to direct
forcefields whose driver works on the arrays of the backend (e.g. the `dummy`
potential), which must then be written with JAX operations. The steps that
need forces from other forcefields run eagerly.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

from functools import partial

from ipi.utils.array_backend import backend
from ipi.utils.cudagraphs import GraphReplay
from ipi.utils.jaxjit import JitReplay
from ipi.utils.jit import JitError
from ipi.utils.messages import verbosity, warning

__all__ = ["replay_steps"]


def unsupported(simulation, system):
    """Returns why the steps of `system` cannot be replayed as CUDA graphs, or
    None if they can."""

    from ipi.engine.barostats import Barostat, BaroBZP, BaroMTK, BaroRGB
    from ipi.engine.motion.dynamics import (
        Dynamics,
        NVEIntegrator,
        NVTIntegrator,
        NPTIntegrator,
        NSTIntegrator,
    )
    from ipi.engine.smotion import Smotion
    from ipi.engine.thermostats import (
        ThermoLangevin,
        ThermoPILE_G,
        ThermoPILE_L,
        ThermoSVR,
        Thermostat,
    )

    motion = system.motion
    if type(motion) is not Dynamics:
        return f"the motion class '{type(motion).__name__}'"
    integrators = (NVEIntegrator, NVTIntegrator, NPTIntegrator, NSTIntegrator)
    if type(motion.integrator) not in integrators:
        return f"the ensemble '{motion.enstype}'"
    thermostats = (Thermostat, ThermoLangevin, ThermoSVR, ThermoPILE_L, ThermoPILE_G)
    if type(motion.thermostat) not in thermostats:
        return f"the thermostat '{type(motion.thermostat).__name__}'"
    if type(motion.barostat) not in (Barostat, BaroBZP, BaroMTK, BaroRGB):
        return f"the barostat '{type(motion.barostat).__name__}'"
    if type(motion.barostat.thermostat) not in (Thermostat, ThermoLangevin):
        name = type(motion.barostat.thermostat).__name__
        return f"the thermostat '{name}' of the barostat"
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
    """Replays the steps of the systems of `simulation` as CUDA graphs, or
    compiles them with JAX, for the systems that support it."""

    for system in simulation.syslist:
        reason = unsupported(simulation, system)
        if reason is None:
            JitStep(system) if backend() == "jax" else GraphStep(system)
        else:
            warning(
                f"The steps are not recorded and replayed: {reason} is not supported",
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
            raise JitError("The random number generator has been replaced")
        super().__call__(step)

    def split(self, bead):
        """Interrupts the recording to compute the forces on `bead` eagerly,
        unless they are already computed, or only the calls to the driver of
        its forcefield are eager."""

        if bead.ufvx_.tainted() and not _eager_driver(bead.ff):
            self.interrupt(partial(_compute, self.components[bead]))


def _compute(component):
    """Computes the forces, potentials and virials of a force component."""

    component.f_.get()
    component.pots_.get()
    component.virs_.get()


class JitStep(JitReplay):
    """Replaces the step of the integrator of a system: runs the first steps
    eagerly, compiles the next one with jax.jit, and calls the compiled step
    at all the following steps."""

    def __init__(self, system):
        integrator = system.motion.integrator
        super().__init__(integrator.step, system.prng)
        integrator.step = self

    def split(self, bead):
        """The forces on `bead` are computed in the compiled step only if the
        driver of its forcefield works on the arrays of the backend."""

        if bead.ufvx_.tainted() and not _eager_driver(bead.ff):
            self.cannot(f"the forces of the forcefield '{bead.ff.name}' are needed")
