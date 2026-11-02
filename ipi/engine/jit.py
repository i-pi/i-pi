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
those it returned while the step was recorded. The calls to drivers that are
`traceable` (e.g. the `dummy` potential) are recorded too, so that the step
is a single CUDA graph.

Only the dynamics that are known to work in this way are replayed (see
`unsupported`), and a `JitError` stops the simulation if a host value read by
a graph changes.

With JAX, the whole step is compiled, including the requests to direct
forcefields whose driver works on the arrays of the backend and can be traced
by jax.jit (`traceable`). The steps that need forces from other forcefields
run eagerly.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

from functools import partial

from ipi.utils.array_backend import backend
from ipi.utils.cudagraphs import GraphReplay
from ipi.utils.jaxjit import JaxReplay
from ipi.utils.jit import JitError
from ipi.utils.messages import verbosity, warning

__all__ = ["compile_steps"]


def unsupported(simulation, system):
    """Returns why the steps of `system` cannot be compiled (replayed as CUDA
    graphs, or compiled with jax.jit), or None if they can."""

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


def compile_steps(simulation):
    """Replays the steps of the systems of `simulation` as CUDA graphs, or
    compiles them with JAX, for the systems that support it."""

    for system in simulation.syslist:
        reason = unsupported(simulation, system)
        if reason is None:
            step = JaxStep if backend() == "jax" else GraphStep
            system.motion.integrator.step = step(system)
        else:
            warning(
                f"The steps are not recorded and replayed: {reason} is not supported",
                verbosity.low,
            )


def _driver_in_step(ff):
    """Whether the requests to `ff` can be made in the compiled step: CUDA
    graphs record them, and call its driver eagerly unless it is `traceable`,
    and jax.jit traces them if the driver is `traceable`."""

    from ipi.engine.forcefields import FFDirect

    return type(ff) is FFDirect and ff.backend_arrays and not ff.threaded


class GraphStep(GraphReplay):
    """The step of the integrator of a system, which runs eagerly for the first
    steps, is recorded as CUDA graphs at the next one, and is replayed for all
    the following steps."""

    def __init__(self, system):
        super().__init__(system.motion.integrator.step, system.prng.torch_generator)
        self.prng = system.prng

    def __call__(self, step=None):
        if (
            self.segments is not None
            and self.prng.torch_generator is not self.generator
        ):
            raise JitError("The random number generator has been replaced")
        super().__call__(step)

    def split(self, bead):
        """Interrupts the recording to compute the forces of the component of
        `bead` eagerly, unless they are already computed (as arrays of the
        backend, from the results of the forcefield, which may be host
        values), or only the calls to the driver of its forcefield are eager."""

        component = bead.component
        computed = (component.f_, component.pots_, component.virs_)
        if any(dep.tainted() for dep in computed) and not _driver_in_step(bead.ff):
            self.interrupt(partial(_compute, component))


def _compute(component):
    """Computes the forces, potentials and virials of a force component."""

    component.f_.get()
    component.pots_.get()
    component.virs_.get()


class JaxStep(JaxReplay):
    """The step of the integrator of a system, which runs eagerly for the first
    steps, is compiled with jax.jit at the next one, and the compiled step runs
    at all the following steps."""

    def __init__(self, system):
        super().__init__(system.motion.integrator.step, system.prng)

    def split(self, bead):
        """The forces on `bead` are computed in the compiled step only if the
        driver of its forcefield works on the arrays of the backend, and can
        be traced by jax.jit."""

        ff = bead.ff
        if bead.component.results_.tainted() and not (
            _driver_in_step(ff) and ff.traceable
        ):
            self.cannot(f"the forces of the forcefield '{ff.name}' are needed")
