"""The behaviour of --jit shared by the CUDA graphs and jax.jit: which dynamics
have their steps compiled (`ipi.engine.jit.unsupported`), and the parameters
fixed in the compiled steps."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import pytest

from ipi_tests.unit_tests.engine.inputs import ISOTROPIC, NO_FORCES, write_input

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

    from ipi.engine.jit import GraphStep, JaxStep

    step = getattr(getattr(sim.syslist[0].motion, "integrator", None), "step", None)
    return isinstance(step, (GraphStep, JaxStep))


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
                "</tau><compressibility units='bar^-1'>4.5e-5</compressibility>"
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


def two_systems(xml, threading, smotion=""):
    """Duplicates the system of an input, with an optional motion of the systems."""

    start, end = xml.index("<system>"), xml.index("</system>") + len("</system>")
    system = xml[start:end]
    systems = system.replace("<system>", "<system prefix='a'>") + system.replace(
        "<system>", "<system prefix='b'>"
    )
    xml = xml[:start] + systems + smotion + xml[end:]
    return xml.replace("threading='false'", f"threading='{threading}'")


@pytest.mark.parametrize(
    "threading,smotion,reason",
    [
        ("true", "", "stepping systems in threads"),
        (
            "false",
            "<smotion mode='remd'><remd><stride>4</stride></remd></smotion>",
            "a motion of multiple systems",
        ),
    ],
)
def test_unsupported_multiple_systems(threading, smotion, reason, simulation):
    """Systems stepped in threads, or with a motion of the systems, are not
    compiled"""

    from ipi.engine.jit import unsupported

    xml = two_systems(write_input("nvt_langevin", NO_FORCES), threading, smotion)
    with simulation(xml) as sim:
        assert unsupported(sim, sim.syslist[0]) == reason


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
