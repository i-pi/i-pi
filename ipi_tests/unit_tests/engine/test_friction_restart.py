"""Continue real XML checkpoints without resetting bath memory or noise."""

import json

import numpy as np
import pytest

from ipi.engine.outputs import CheckpointOutput
from ipi.scripting import InteractiveSimulation
from ipi.utils.io import NumpyEncoder

DRIVER = """
import json
import numpy as np
from ipi.pes.dummy import Dummy_driver
__DRIVER_CLASS__ = "RestartDriver"
class RestartDriver(Dummy_driver):
    def compute_structure(self, cell, pos):
        x, y, z = pos.reshape(-1)
        j = np.array([[1+x*x, 0., 0.], [0., 1., 0.]])
        extra = {"friction_coupling": [x+x**3/3, y],
                 "friction_coupling_jacobian": j.tolist(),
                 "friction_gamma": (j.T@j).tolist(),
                 "friction_meta": {"schema": "ipi_friction_v1",
                    "active_atoms": [0], "channel_labels": ["x", "y"],
                    "jacobian_units": "ipi_internal_mass_scaled",
                    "gamma_units": "ipi_internal_mass_scaled"}}
        return .5*np.sum(pos**2), -pos, np.zeros((3,3)), json.dumps(extra)
"""


def advance(sim):
    sim.get_motion_step()(step=sim.step)
    sim.step += 1


def state(sim):
    system = sim.syslist[0]
    friction = system.motion.friction
    # A restored variable bath defers allocation until the driver is queried.
    if friction.variable_friction and friction.bath.s is None:
        friction.bath._get_non_markovian_jacobian()
    return [
        system.beads.q.copy(),
        system.beads.p.copy(),
        friction.bath.s.copy(),
        np.array(friction.ediss),
    ]


def make_simulation(tmp_path, variable, splitting, mode="nve-f", doublewell=False):
    driver = tmp_path / "restart_driver.py"
    driver_code = DRIVER
    if doublewell:
        driver_code = driver_code.replace(
            ".5*np.sum(pos**2), -pos", "np.sum(.25*pos**4-.5*pos**2), pos-pos**3"
        )
    driver.write_text(driver_code)
    thermostat = (
        '<thermostat mode="langevin"><tau>1.</tau></thermostat>'
        if mode == "nvt-f"
        else ""
    )
    matrices = np.array(
        [[[0.0, c], [-c, a]] for c, a in [(0.2, 0.3), (0.3, 0.4), (0.4, 0.5)]]
    )
    ap = ", ".join(map(str, matrices.ravel()))
    xml = f"""<simulation verbosity="quiet" floatformat="%.17e">
    <prng><seed>9876</seed></prng><total_steps>1000</total_steps>
    <ffdirect name="test" pbc="false"><pes>custom</pes><pes_path>{driver}</pes_path></ffdirect>
    <output prefix="{tmp_path / 'out'}"></output>
    <system><beads natoms="1" nbeads="3">
      <q shape="(3,3)">[.1,.2,.3,.3,.2,.1,.2,.3,.1]</q>
      <p shape="(3,3)">[.2,.1,.3,.1,.3,.2,.3,.2,.1]</p>
      <m shape="(1)">[100]</m><names shape="(1)">[H]</names>
    </beads><cell shape="(3,3)">[10,0,0,0,10,0,0,0,10]</cell>
    <ensemble><temperature>0.01</temperature></ensemble>
    <forces><force forcefield="test"/></forces>
    <motion mode="dynamics"><fixcom>False</fixcom><dynamics mode="{mode}" splitting="{splitting}">{thermostat}
      <timestep>0.05</timestep><friction>
        <bath_mode>non-markovian</bath_mode><variable_friction>{variable}</variable_friction>
        <Ap shape="(3,2,2)">[{ap}]</Ap>
        <Lambda shape="(4,2)">[0,.5,1,.3,2,.1,10,0]</Lambda>
      </friction></dynamics></motion></system></simulation>"""
    return InteractiveSimulation(xml)


@pytest.mark.parametrize("splitting", ["fbabf", "symmetric-gle"])
@pytest.mark.parametrize("variable", [False, True])
@pytest.mark.parametrize("split_step", [0, 7])
def test_checkpoint_continues_same_stochastic_trajectory(
    tmp_path, variable, split_step, splitting
):
    original = make_simulation(tmp_path, variable, splitting)
    for _ in range(split_step):
        advance(original)
    path = tmp_path / "RESTART"
    checkpoint = CheckpointOutput(filename=str(path), stride=1)
    checkpoint.bind(original)
    # Interactive stepping already advanced the next-step counter.
    checkpoint.write(store=False)
    saved_rng = json.dumps(original.prng.state, cls=NumpyEncoder)
    expected = []
    for _ in range(10):
        advance(original)
        expected.append(state(original))
    final_rng = json.dumps(original.prng.state, cls=NumpyEncoder)

    resumed = InteractiveSimulation(path.read_text())
    assert resumed.step == split_step
    assert json.dumps(resumed.prng.state, cls=NumpyEncoder) == saved_rng
    for reference in expected:
        advance(resumed)
        for actual, wanted in zip(state(resumed), reference):
            np.testing.assert_allclose(actual, wanted, rtol=1e-12, atol=1e-14)
    assert json.dumps(resumed.prng.state, cls=NumpyEncoder) == final_rng


def test_restart_rejects_invalid_auxiliary_state():
    from ipi.engine.friction import Friction

    with pytest.raises(ValueError, match="finite three-dimensional"):
        Friction(s=np.full((2, 1, 1), np.nan))
    with pytest.raises(ValueError, match="non-markovian"):
        Friction(bath_mode="markovian", s=np.ones((1, 1, 1)))


@pytest.mark.parametrize("variable", [False, True])
@pytest.mark.parametrize("splitting", ["symmetric-gle"])
@pytest.mark.parametrize("doublewell", [False, True])
def test_gle_matches_independent_update(tmp_path, variable, splitting, doublewell):
    """Numerical oracle using direct bead/mode algebra, not integrator methods."""
    from scipy.linalg import expm
    from ipi.utils.nmtransform import mk_nm_matrix

    sim = make_simulation(tmp_path, variable, splitting, doublewell=doublewell)
    system = sim.syslist[0]
    friction = system.motion.friction
    bath = friction.bath
    if variable:
        bath._get_non_markovian_jacobian()
    bath.s[:] = np.arange(bath.s.size).reshape(bath.s.shape) * 0.02 + 0.1
    sim.prng.gvec = lambda shape: np.zeros(shape)  # Deterministic OU drift.
    q, p, s, _ = state(sim)
    h = float(system.motion.dt)
    mass = np.asarray(system.beads.m3).copy()
    c = mk_nm_matrix(3)
    alpha = friction.alpha.copy()
    theta = bath.theta.copy()
    ou = np.array([expm(-a * h / 2) for a in bath.A_aux])

    def kick(q, p, s):
        if variable:
            f = np.c_[q[:, 0] + q[:, 0] ** 3 / 3, q[:, 1]]
            jac = np.zeros((3, 2, 3))
            jac[:, 0, 0], jac[:, 1, 1] = 1 + q[:, 0] ** 2, 1.0
        else:
            f = q.copy()
            jac = np.broadcast_to(np.eye(3), (3, 3, 3))
        mf = -np.einsum("bci,bc->bi", jac, c.T @ (alpha[:, None] * (c @ f)))
        bath_force = -np.einsum(
            "bci,bc->bi", jac, c.T @ np.einsum("ma,mac->mc", theta, s)
        )
        physical_force = q - q**3 if doublewell else -q
        return p + h / 2 * (physical_force + mf + bath_force), jac

    # O/2, B_P,V/2 and B_P,F/2.
    s = np.einsum("mab,mbc->mac", ou, s)
    p, jac = kick(q, p, s)
    if splitting == "symmetric-gle":
        drive = c @ np.einsum("bci,bi->bc", jac, p / mass)
        s += h / 2 * theta[:, :, None] * drive[:, None, :]
    # A_omega(h): analytic full-step free-ring-polymer propagation.
    qn, pn = c @ q, c @ p
    w = np.asarray(system.nm.omegak)
    cos = np.cos(w * h)[:, None]
    sin_over_w = (h * np.sinc(w * h / np.pi))[:, None]
    q = c.T @ (cos * qn + sin_over_w * pn / mass)
    p = c.T @ (cos * pn - (w * w)[:, None] * sin_over_w * mass * qn)
    # B_s(h), using the NEW coordinates and velocities.
    if variable:
        jac = np.zeros((3, 2, 3))
        jac[:, 0, 0], jac[:, 1, 1] = 1 + q[:, 0] ** 2, 1.0
    else:
        jac = np.broadcast_to(np.eye(3), (3, 3, 3))
    drive = c @ np.einsum("bci,bi->bc", jac, p / mass)
    drive_dt = h / 2 if splitting == "symmetric-gle" else h
    s += drive_dt * theta[:, :, None] * drive[:, None, :]
    p, _ = kick(q, p, s)
    s = np.einsum("mab,mbc->mac", ou, s)

    advance(sim)
    for actual, expected in zip(state(sim)[:3], [q, p, s]):
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)


def test_restart_rejects_auxiliary_shape_mismatch():
    from types import SimpleNamespace
    from ipi.engine.friction import Friction

    friction = Friction(
        variable_friction=False,
        s=np.ones((1, 1, 2)),
        Ap=np.array([[0.0, 0.2], [-0.2, 0.3]]),
    )
    friction.beads = SimpleNamespace(nbeads=1, natoms=1)
    with pytest.raises(ValueError, match="restart s has shape"):
        friction._ensure_bath_bound()


def test_pending_restart_state_survives_input_roundtrip():
    from ipi.engine.friction import Friction
    from ipi.inputs.friction import InputFriction

    signature = [
        [0, 1, 2],
        ["x", "y"],
        "ipi_internal_mass_scaled",
        "ipi_internal_mass_scaled",
    ]
    source = Friction(
        s=np.arange(6.0).reshape(3, 1, 2),
        ediss=1.25,
        bath_contract=json.dumps(signature),
    )
    stored = InputFriction()
    stored.store(source)  # No driver query or random draw needed to store it.
    source._initial_s[:] = -10  # Checkpoint arrays must own their data.
    restored = stored.fetch()
    np.testing.assert_array_equal(restored._initial_s, np.arange(6.0).reshape(3, 1, 2))
    assert restored.ediss == 1.25
    assert restored._friction_contract_signature == (
        tuple(signature[0]),
        tuple(signature[1]),
        *signature[2:],
    )


@pytest.mark.parametrize(
    "field,value,message",
    [
        ("bath_mode", "markovian", "non-markovian"),
        ("nmts", [2], "single force timestep"),
        ("mode", "manual", "normal_modes"),
        ("propagator", "cayley", "normal_modes"),
        ("open_paths", [0], "closed"),
        ("bosons", [0], "distinguishable"),
        ("fixcom", True, "fixcom=False"),
    ],
)
def test_symmetric_gle_rejects_unsupported_settings(field, value, message):
    from types import SimpleNamespace
    from ipi.engine.motion.dynamics import NVEIntegratorWithFriction

    integrator = SimpleNamespace(
        friction=SimpleNamespace(bath_mode="non-markovian"),
        nmts=[1],
        nm=SimpleNamespace(mode="rpmd", propagator="exact", open_paths=[], bosons=[]),
        fixcom=False,
        fixatoms_dof=[],
    )
    target = integrator
    if field == "bath_mode":
        target = integrator.friction
    elif field in ("mode", "propagator", "open_paths", "bosons"):
        target = integrator.nm
    setattr(target, field, value)
    with pytest.raises(ValueError, match=message):
        NVEIntegratorWithFriction._validate_symmetric_gle(integrator)
