"""The behaviour of --jit shared by the CUDA graphs and jax.jit: which dynamics
have their steps compiled (`ipi.engine.jit.unsupported`), and the parameters
fixed in the compiled steps."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi_tests.unit_tests.engine.inputs import (
    CASES,
    DUMMY,
    ISOTROPIC,
    NO_FORCES,
    PROPERTIES,
    TOLERANCE,
    write_input,
)

try:
    import torch

    cuda_available = torch.cuda.is_available()
except ImportError:
    cuda_available = False

# the options of the compiled simulations, with either backend
JIT = [
    pytest.param(dict(backend="jax", jit=True), id="jax"),
    pytest.param(
        dict(backend="torch", device="cuda", jit=True),
        id="cuda",
        marks=pytest.mark.skipif(not cuda_available, reason="CUDA not available"),
    ),
]


def compiled(sim):
    """Whether the steps of the first system are compiled."""

    from ipi.engine.jit import compiled_step

    return compiled_step(sim.syslist[0]) is not None


@pytest.mark.parametrize(
    "case,change,reason",
    [
        ("nve", None, None),
        ("nvt_svr", None, None),
        ("pimd_pile_g", None, None),
        ("nvt_gle", None, "the thermostat 'ThermoGLE'"),
        ("npt_flexible", None, None),
        (
            "npt_isotropic",
            (ISOTROPIC, ISOTROPIC.replace("'langevin'", "'svr'")),
            "the thermostat 'ThermoSVR' of the barostat",
        ),
        (
            "npt_isotropic",
            (
                ISOTROPIC,
                "<barostat mode='stochastic-rescaling'><tau units='femtosecond'>200"
                "</tau><bulk_modulus units='gigapascal'>2.2</bulk_modulus>"
                "</barostat>",
            ),
            "the barostat 'BaroSCR'",
        ),
        ("nvt_langevin", ("mode='dynamics'", "mode='dummy'"), "the motion class"),
        ("pimd_nvt", ("mode='nvt'", "mode='nvt-cc'"), "the ensemble 'nvt-cc'"),
        (
            "nvt_langevin",
            (
                "<motion mode='dynamics'>",
                "<motion mode='dynamics'><fixatoms>[0]</fixatoms>",
            ),
            "fixed atoms",
        ),
        (
            "pimd_nvt",
            (
                "<motion mode='dynamics'>",
                "<normal_modes><open_paths>[0]</open_paths></normal_modes>"
                "<motion mode='dynamics'>",
            ),
            "open paths or bosons",
        ),
    ],
)
def test_unsupported_steps(case, change, reason, simulation):
    """The steps of the supported dynamics are compiled, and the reason why the
    others are not is given"""

    from ipi.engine.jit import compile_steps, unsupported

    xml = write_input(case, NO_FORCES)
    if change is not None:
        xml = xml.replace(*change)
    with simulation(xml) as sim:
        if reason is None:
            assert unsupported(sim, sim.syslist[0]) is None
        else:
            assert reason in unsupported(sim, sim.syslist[0])
            compile_steps(sim)
            assert not compiled(sim)


# exchanges of the replicas tried at every step
REMD = "<smotion mode='remd'><remd><stride>1</stride></remd></smotion>"


def replicas(xml, temperatures, smotion=REMD, interleave="true"):
    """Replicates the system of an input at the given temperatures, with an
    optional motion of the systems."""

    start, end = xml.index("<system>"), xml.index("</system>") + len("</system>")
    system = xml[start:end]
    systems = "".join(
        system.replace("<system>", f"<system prefix='r{i}'>").replace(
            "<temperature units='kelvin'>300</temperature>",
            f"<temperature units='kelvin'>{temperature}</temperature>",
        )
        for i, temperature in enumerate(temperatures)
    )
    xml = xml[:start] + systems + smotion + xml[end:]
    return xml.replace("interleave='false'", f"interleave='{interleave}'")


@pytest.mark.parametrize(
    "case,smotion,change,reason",
    [
        ("nvt_langevin", "", None, None),
        ("pimd_nvt", REMD, None, None),
        (
            "nvt_langevin",
            "<smotion mode='multi'>" + REMD + "</smotion>",
            None,
            "a motion of multiple systems",
        ),
        ("npt_isotropic", REMD, None, "replica exchange with a barostat"),
        (
            "nvt_langevin",
            REMD,
            (
                "<temperature units='kelvin'>310</temperature>",
                "<temperature "
                "units='kelvin'>310</temperature><hamiltonian_weights>[0.5]"
                "</hamiltonian_weights>",
            ),
            "replica exchange of the weights of the forces",
        ),
    ],
)
def test_supported_multiple_systems(case, smotion, change, reason, simulation):
    """The steps of several systems, which are stepped in turn, are compiled,
    also with replica exchange, unless it exchanges parameters that are fixed
    in the compiled steps"""

    from ipi.engine.jit import unsupported

    xml = replicas(write_input(case, NO_FORCES), [300, 310], smotion)
    if change is not None:
        xml = xml.replace(*change)
    with simulation(xml) as sim:
        assert unsupported(sim, sim.syslist[0]) == reason


def record_batches(forcefield):
    """Wraps the driver of a direct forcefield, and returns the list of the
    numbers of structures that it is given at each call."""

    batches = []
    driver = forcefield.driver

    def call(cell, pos):
        batches.append(len(cell) if isinstance(cell, list) else 1)
        return driver(cell, pos)

    forcefield.driver = call
    return batches


@pytest.mark.parametrize("options", JIT)
@pytest.mark.parametrize("batched", [False, True])
@pytest.mark.parametrize("case", ["nvt_langevin", "pimd_nvt"])
def test_replica_exchange_is_compiled(
    case, batched, options, simulation, device_random
):
    """The steps of replicas whose temperatures are exchanged at every step,
    also before the steps are compiled, are compiled together, and follow the
    eager ones. A batched driver evaluates the structures of all the replicas
    at once, also in the compiled steps, which repeat the recorded calls"""

    from ipi.engine.jit import compiled_step

    pytest.importorskip(options["backend"])
    temperatures = [300, 330, 360, 400]
    nbatch = len(temperatures) * CASES[case][0] if batched else 1
    forcefield = (
        f"<ffdirect name='ff'><pes>dummy</pes><batch_size>{nbatch}</batch_size></ffdirect>",
        "<force forcefield='ff'/>",
    )
    xml = replicas(write_input(case, forcefield), temperatures)
    eager_options = {k: v for k, v in options.items() if k != "jit"}
    with simulation(xml, **eager_options) as eager, simulation(xml, **options) as sim:
        batches = record_batches(sim.fflist["ff"])
        temperatures = set()
        for step in range(10):
            eager.run_step(step)
            sim.run_step(step)
            temperatures.add(float(sim.syslist[0].ensemble.temp))
        assert all(compiled_system(s) for s in sim.syslist)
        assert len({compiled_step(s) for s in sim.syslist}) == 1
        assert set(batches) == {nbatch}
        assert list(sim.smotion.repindex) == list(eager.smotion.repindex)
        # the replica has been exchanged
        assert len(temperatures) > 1
        for name in PROPERTIES:
            np.testing.assert_allclose(
                sim.properties(name), eager.properties(name), **TOLERANCE
            )
        for a, b in zip(eager.syslist, sim.syslist):
            for x, y in [(a.beads.q, b.beads.q), (a.beads.p, b.beads.p)]:
                np.testing.assert_allclose(ab.to_numpy(y), ab.to_numpy(x), **TOLERANCE)


@pytest.mark.parametrize("options", JIT)
@pytest.mark.parametrize(
    "interleave,forcefield,together",
    [("true", DUMMY, True), ("false", DUMMY, False), ("true", NO_FORCES, False)],
)
def test_steps_compiled_together(interleave, forcefield, together, options, simulation):
    """The steps of several systems are compiled together if the systems are
    stepped as tasks, and all their forces are computed in the step"""

    from ipi.engine.jit import compiled_step

    pytest.importorskip(options["backend"])
    xml = replicas(write_input("nvt_langevin", forcefield), [300, 310], "", interleave)
    with simulation(xml, **options) as sim:
        steps = [compiled_step(s) for s in sim.syslist]
        assert all(step is not None for step in steps)
        assert (steps[0] is steps[1]) == together


@pytest.mark.parametrize("options", JIT)
@pytest.mark.parametrize("remd", [False, True])
def test_stale_objects_are_refreshed_only_with_exchanges(
    remd, options, simulation, device_random, monkeypatch
):
    """Before the eager steps and the compilation, the objects that the steps
    read are brought up to date if they are stale, which happens only when
    something else modifies the parameters of the steps, e.g. the exchanges
    of the replicas: the compiled steps of other simulations do not change"""

    from ipi.utils.jit import Replay

    pytest.importorskip(options["backend"])
    stale = []
    refresh = Replay.refresh

    def counting_refresh(self):
        stale.append(sum(dep.tainted() for dep in self.computed.values()))
        refresh(self)

    monkeypatch.setattr(Replay, "refresh", counting_refresh)
    smotion = REMD if remd else ""
    xml = replicas(write_input("pimd_nvt", DUMMY), [300, 330], smotion)
    with simulation(xml, **options) as sim:
        for step in range(8):
            sim.run_step(step)
        assert all(compiled_system(s) for s in sim.syslist)
    assert stale
    assert any(stale) == remd


def compiled_system(system):
    """Whether the steps of a system are compiled."""

    from ipi.engine.jit import compiled_step

    step = compiled_step(system)
    return step is not None and step.compiled is not None


@pytest.mark.parametrize("options", JIT)
@pytest.mark.parametrize(
    "case,owner,name",
    [
        ("nvt_langevin", "motion", "dt"),
        ("npt_isotropic", "ensemble", "pext"),
        ("npt_isotropic", "ensemble", "temp"),
    ],
)
def test_changed_parameters_stop_the_compiled_steps(
    case, owner, name, options, simulation, device_random
):
    """The time step, and the pressure and temperature of the barostats, are
    fixed in the compiled steps"""

    from ipi.utils.jit import JitError

    pytest.importorskip(options["backend"])
    with simulation(write_input(case, NO_FORCES), **options) as sim:
        for step in range(4):
            sim.run_step(step)
        assert compiled(sim)
        obj = getattr(sim.syslist[0], owner)
        setattr(obj, name, getattr(obj, name) * 1.1)
        with pytest.raises(JitError, match="cannot change"):
            sim.run_step(4)
