"""The petjax driver, with neighbor lists from vesin or computed with JAX, on a
small randomly initialised PET model."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("petjax")

import jax  # noqa: E402
from ase import Atoms  # noqa: E402
from ase.io import write  # noqa: E402

from ipi.pes.petjax import PETJAX_driver  # noqa: E402
from ipi.utils.units import unit_to_internal  # noqa: E402

# a tiny model, with an adaptive cutoff and the direct force and stress heads
CONFIG = dict(
    d_pet=8,
    d_node=16,
    d_head=8,
    d_feedforward=8,
    num_heads=2,
    num_attention_layers=1,
    num_gnn_layers=2,
    cutoff=3.0,
    cutoff_width=1.0,
    cutoff_width_adaptive=1.0,
    adaptive_cutoff_method="solver",
    num_neighbors_adaptive=4,
    max_atomic_number=8,
)
SHIFTS = {1: -0.5, 8: -2.0}
BOHR = unit_to_internal("length", "angstrom", 1.0)


@pytest.fixture(scope="module")
def checkpoint(tmp_path_factory):
    """A pet-jax checkpoint of the tiny model, with random weights."""

    import jax.numpy as jnp
    from marathon.io import write_msgpack, write_yaml
    from petjax import UPET

    positions = jnp.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.2, 0.0]])
    pairs = [(i, j) for i in range(3) for j in range(3) if i != j]
    centers = jnp.array([i for i, _ in pairs])
    neighbors = jnp.array([j for _, j in pairs])
    inputs = dict(
        R_ij=positions[neighbors] - positions[centers],
        centers=centers,
        neighbors=neighbors,
        atomic_numbers=jnp.array([8, 1, 1]),
        reverse=jnp.array([pairs.index((j, i)) for i, j in pairs]),
        pair_mask=jnp.ones(6, dtype=bool),
        atom_mask=jnp.ones(3, dtype=bool),
    )
    model = UPET(**CONFIG, direct_forces=True, direct_stress=True)
    variables = model.init(jax.random.key(0), **inputs)
    path = tmp_path_factory.mktemp("petjax")
    write_msgpack(path / "model.msgpack", variables)
    write_yaml(path / "metadata.yaml", {"config": CONFIG, "shifts": SHIFTS})
    return str(path)


def water(n, cell, pbc, seed=0):
    """`n` water-like molecules at random positions in `cell`."""

    rng = np.random.default_rng(seed)
    centers = rng.random((n, 3)) @ np.asarray(cell)
    offsets = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])
    positions = (centers[:, None, :] + offsets).reshape(-1, 3)
    return Atoms("OH2" * n, positions=positions, cell=cell, pbc=pbc)


STRUCTURES = {
    # a triclinic cell shorter than the cutoff along one axis, so that the
    # pairs come from several cell images
    "periodic": water(4, [[4.5, 0.0, 0.0], [0.8, 4.2, 0.0], [0.3, -0.4, 2.6]], True),
    "molecule": water(2, np.eye(3) * 12.0, False),
}


def template(tmp_path, name):
    path = str(tmp_path / f"{name}.xyz")
    write(path, STRUCTURES[name])
    return path


def ipi_inputs(atoms, displacement=0.0, seed=1):
    """The cell (vectors as columns) and positions of `atoms`, in atomic units,
    with the positions displaced at random and possibly outside the cell."""

    rng = np.random.default_rng(seed)
    positions = atoms.positions + rng.normal(scale=displacement, size=(len(atoms), 3))
    return np.array(atoms.cell).T * BOHR, positions * BOHR


def as_numpy(result):
    energy, forces, virial, _ = result
    return (
        float(energy),
        np.asarray(forces, dtype=float).reshape(-1, 3),
        np.asarray(virial, dtype=float).reshape(3, 3),
    )


@pytest.mark.parametrize("non_conservative", [False, True])
@pytest.mark.parametrize("name", list(STRUCTURES))
def test_jax_neighbors_match_vesin(name, non_conservative, checkpoint, tmp_path):
    """The neighbor lists computed with JAX give the same energy, forces and
    virial as those of vesin, also for atoms outside the cell"""

    xyz = template(tmp_path, name)
    options = dict(non_conservative=non_conservative)
    vesin = PETJAX_driver(xyz, checkpoint, neighbors="vesin", **options)
    jax_lists = PETJAX_driver(xyz, checkpoint, neighbors="jax", **options)
    atoms = STRUCTURES[name].copy()
    # moves one molecule by a cell vector
    atoms.positions[:3] += np.array(atoms.cell)[0] * atoms.pbc[0]
    cell, positions = ipi_inputs(atoms, displacement=0.05)
    e1, f1, v1 = as_numpy(vesin.compute(cell, positions))
    e2, f2, v2 = as_numpy(jax_lists.compute(cell, positions))
    # float32 model, summed in different orders
    assert e2 == pytest.approx(e1, abs=1e-6)
    np.testing.assert_allclose(f2, f1, rtol=2e-5, atol=1e-6)
    np.testing.assert_allclose(v2, v1, rtol=2e-5, atol=1e-6)


def test_batches_match_single_structures(checkpoint, tmp_path):
    """A batch is evaluated as its structures one at a time"""

    driver = PETJAX_driver(template(tmp_path, "periodic"), checkpoint, neighbors="jax")
    inputs = [ipi_inputs(STRUCTURES["periodic"], 0.05, seed) for seed in range(3)]
    batch = driver.compute([c for c, _ in inputs], [p for _, p in inputs])
    for (cell, positions), result in zip(inputs, batch):
        single = as_numpy(driver.compute(cell, positions))
        for a, b in zip(as_numpy(result), single):
            np.testing.assert_allclose(a, b, rtol=1e-5, atol=1e-6)


def test_forces_and_virial_are_derivatives_of_the_energy(checkpoint, tmp_path):
    """The forces and the virial are minus the derivatives of the energy with
    respect to the positions and to a strain of the cell, in atomic units"""

    driver = PETJAX_driver(
        template(tmp_path, "periodic"), checkpoint, neighbors="jax", dtype="float64"
    )
    cell, positions = ipi_inputs(STRUCTURES["periodic"], displacement=0.05)
    _, forces, virial = as_numpy(driver.compute(cell, positions))
    h = 1e-5

    def energy(c, p):
        return as_numpy(driver.compute(c, p))[0]

    for atom, axis in [(0, 0), (4, 2), (7, 1)]:
        dp = np.zeros_like(positions)
        dp[atom, axis] = h
        derivative = (energy(cell, positions + dp) - energy(cell, positions - dp)) / (
            2 * h
        )
        assert forces[atom, axis] == pytest.approx(-derivative, rel=1e-5, abs=1e-8)
    for a, b in [(0, 0), (0, 1), (2, 1)]:
        strain = np.zeros((3, 3))
        strain[a, b] = strain[b, a] = h
        plus, minus = np.eye(3) + strain, np.eye(3) - strain
        derivative = (
            energy(plus @ cell, positions @ plus.T)
            - energy(minus @ cell, positions @ minus.T)
        ) / (2 * h)
        expected = derivative if a == b else derivative / 2
        assert virial[a, b] == pytest.approx(-expected, rel=1e-5, abs=1e-8)


def test_exceeded_capacity_stops_the_evaluation(checkpoint, tmp_path):
    """Neighbor lists computed with JAX that cannot hold all the pairs within
    the cutoff stop the evaluation with an error"""

    driver = PETJAX_driver(
        template(tmp_path, "periodic"), checkpoint, neighbors="jax", pair_margin=0.5
    )
    cell, positions = ipi_inputs(STRUCTURES["periodic"])
    with pytest.raises(Exception, match="increase pair_margin"):
        jax.block_until_ready(driver.compute(cell, positions))


def test_steps_are_compiled_with_jax(checkpoint, tmp_path, simulation):
    """With the JAX backend and --jit, the steps of a PIMD simulation, including
    the evaluation of the model on all the beads, are compiled into one call,
    and follow the eager steps"""

    from ipi.engine.jit import compiled_step
    from ipi.utils.array_backend import to_numpy

    pytest.importorskip("jax")
    xyz = template(tmp_path, "periodic")
    xml = f"""<simulation verbosity='quiet' safe_stride='100'>
<output prefix='x'/><total_steps>100</total_steps><prng><seed>3</seed></prng>
<ffdirect name='pet' pbc='false'><pes>petjax</pes>
<parameters>{{template:{xyz},checkpoint:{checkpoint},neighbors:jax}}</parameters>
<batch_size>2</batch_size></ffdirect>
<system><initialize nbeads='2'><file mode='ase' units='angstrom'>{xyz}</file>
<velocities mode='thermal' units='kelvin'>300</velocities></initialize>
<forces><force forcefield='pet'/></forces>
<motion mode='dynamics'><dynamics mode='nvt'><timestep units='femtosecond'>0.5</timestep>
<thermostat mode='pile_l'><tau units='femtosecond'>100</tau></thermostat></dynamics></motion>
<ensemble><temperature units='kelvin'>300</temperature></ensemble></system></simulation>"""
    with simulation(xml, backend="jax") as eager, simulation(
        xml, backend="jax", jit=True
    ) as compiled:
        for step in range(5):
            eager.run_step(step)
            compiled.run_step(step)
        assert compiled_step(compiled.syslist[0]).compiled is not None
        a, b = eager.syslist[0].beads, compiled.syslist[0].beads
        np.testing.assert_allclose(to_numpy(b.q), to_numpy(a.q), rtol=0, atol=1e-8)
        np.testing.assert_allclose(to_numpy(b.p), to_numpy(a.p), rtol=0, atol=1e-8)


def test_pet_mad_matches_metatomic(request, tmp_path):
    """PET-MAD-xs through pet-jax and through metatomic gives the same energy,
    forces and virial, to the precision of float32"""

    checkpoint = request.config.getoption("--ipi-petjax-checkpoint")
    model = request.config.getoption("--ipi-metatomic-model")
    if not checkpoint or not model:
        pytest.skip("needs --ipi-petjax-checkpoint and --ipi-metatomic-model")
    from ase.io import read

    from ipi.pes.metatomic import MetatomicDriver

    xyz = str(
        Path(__file__).parents[3] / "examples" / "clients" / "petjax" / "init.xyz"
    )
    atoms = read(xyz)
    cell, positions = ipi_inputs(atoms, displacement=0.02)
    reference = as_numpy(MetatomicDriver(xyz, model).compute(cell, positions))
    for neighbors in ("vesin", "jax"):
        driver = PETJAX_driver(xyz, checkpoint, neighbors=neighbors)
        energy, forces, virial = as_numpy(driver.compute(cell, positions))
        assert energy == pytest.approx(reference[0], abs=1e-4)
        np.testing.assert_allclose(forces, reference[1], atol=1e-4)
        np.testing.assert_allclose(virial, reference[2], atol=1e-3)
