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
a graph changes. Replica exchange is supported as long as the parameters that
it exchanges reach the steps as arrays: the temperature of the thermostats and
of the normal modes, but not that of the barostats, their pressure, or the
weights of the forces.

With `interleave`, the steps of several systems whose forces are all computed
in the step are compiled together (see `_Together`): their integrator steps
run as tasks while they are recorded, so that a direct forcefield evaluates
the requests of all the systems in the same batches, which the compiled step
keeps. Otherwise the steps of each system are compiled on their own, and the
systems are stepped in turn.

With JAX, the whole step is compiled, including the requests to direct
forcefields whose driver works on the arrays of the backend and can be traced
by jax.jit (`traceable`). The steps that need forces from other forcefields
run eagerly.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

from functools import partial

import numpy as np

from ipi.utils import scheduler
from ipi.utils.array_backend import backend
from ipi.utils.cudagraphs import GraphReplay
from ipi.utils.jaxjit import JaxReplay
from ipi.utils.jit import JitError
from ipi.utils.messages import verbosity, warning

__all__ = ["compile_steps", "compiled_step"]


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
    from ipi.engine.smotion import Smotion, ReplicaExchange
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
    smotions = (Smotion, ReplicaExchange)
    if simulation.smotion is not None and type(simulation.smotion) not in smotions:
        return "a motion of multiple systems"
    if type(simulation.smotion) is ReplicaExchange:
        # the barostats and the sums of the forces read the parameters that
        # would be exchanged as host values, fixed in the compiled steps
        if type(motion.barostat) is not Barostat:
            return "replica exchange with a barostat"
        weights = {
            (tuple(np.ravel(s.ensemble.hweights)), tuple(np.ravel(s.ensemble.bweights)))
            for s in simulation.syslist
        }
        if len(weights) > 1:
            return "replica exchange of the weights of the forces"
    return None


def compile_steps(simulation):
    """Replays the steps of the systems of `simulation` as CUDA graphs, or
    compiles them with JAX, for the systems that support it: together, if
    they are several, stepped as tasks (`interleave`), and all their forces are
    computed in the step."""

    systems = []
    for system in simulation.syslist:
        reason = unsupported(simulation, system)
        if reason is None:
            systems.append(system)
        else:
            warning(
                f"The steps are not recorded and replayed: {reason} is not supported",
                verbosity.low,
            )
    step = JaxStep if backend() == "jax" else GraphStep
    if (
        len(systems) > 1
        and len(systems) == len(simulation.syslist)
        and simulation.interleave
        and all(_driver_in_step(c.ff) for s in systems for c in s.forces.mforces)
    ):
        together = _Together(step(systems), len(systems))
        for index, system in enumerate(systems):
            system.motion.integrator.step = partial(together.step, index)
    else:
        for system in systems:
            system.motion.integrator.step = step([system])


def compiled_step(system):
    """Returns the compiled step (a GraphStep or a JaxStep) that makes the
    step of `system`, and possibly those of other systems, or None."""

    step = getattr(getattr(system.motion, "integrator", None), "step", None)
    together = getattr(getattr(step, "func", None), "__self__", None)
    if isinstance(together, _Together):
        return together.replay
    return step if isinstance(step, (GraphStep, JaxStep)) else None


def _integrators(systems):
    """Returns the function that makes the integrator steps of `systems`: those
    of several systems run as tasks, which switch while they wait for forces,
    so that the requests of all of them are evaluated together."""

    steps = [system.motion.integrator.step for system in systems]
    if len(steps) == 1:
        return steps[0]

    def step(step=None):
        scheduler.run([partial(f, step) for f in steps])

    return step


class _Together:
    """The integrator steps of several systems, compiled together. The first of
    them that is called at a step of the simulation makes the steps of all the
    systems, and the others then return, as their step has been made.

    Attributes:
        replay: The compiled step of all the systems.
        calls: The number of calls of the step of each system.
        done: The number of steps made by the compiled step.
    """

    def __init__(self, replay, nsystems):
        self.replay = replay
        self.calls = [0] * nsystems
        self.done = 0

    def step(self, index, step=None):
        """The integrator step of the system `index`."""

        self.calls[index] += 1
        if self.calls[index] > self.done:
            self.done += 1
            self.replay(step)


def _driver_in_step(ff):
    """Whether the requests to `ff` can be made in the compiled step: CUDA
    graphs record them, and call its driver eagerly unless it is `traceable`,
    and jax.jit traces them if the driver is `traceable`."""

    from ipi.engine.forcefields import FFDirect

    return type(ff) is FFDirect and ff.backend_arrays


class GraphStep(GraphReplay):
    """The step of the integrators of some systems, which runs eagerly for the
    first steps, is recorded as CUDA graphs at the next one, and is replayed
    for all the following steps. The systems share the random number
    generator of the simulation."""

    def __init__(self, systems):
        super().__init__(_integrators(systems), systems[0].prng.torch_generator)
        self.prng = systems[0].prng
        self.together = len(systems) > 1

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
            if self.together:
                # the other systems would run, unrecorded, while it waits
                raise JitError(
                    f"The forces of the forcefield '{bead.ff.name}' cannot be "
                    "computed eagerly in the steps of several systems"
                )
            self.interrupt(partial(_compute, component))


def _compute(component):
    """Computes the forces, potentials and virials of a force component."""

    component.f_.get()
    component.pots_.get()
    component.virs_.get()


class JaxStep(JaxReplay):
    """The step of the integrators of some systems, which runs eagerly for the
    first steps, is compiled with jax.jit at the next one, and the compiled
    step runs at all the following steps. The systems share the random number
    generator of the simulation."""

    def __init__(self, systems):
        super().__init__(_integrators(systems), systems[0].prng)

    def split(self, bead):
        """The forces on `bead` are computed in the compiled step only if the
        driver of its forcefield works on the arrays of the backend, and can
        be traced by jax.jit."""

        ff = bead.ff
        if bead.component.results_.tainted() and not (
            _driver_in_step(ff) and ff.traceable
        ):
            self.cannot(f"the forces of the forcefield '{ff.name}' are needed")
