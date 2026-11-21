"""Combined ring-forcefield contract, using the ordinary force-request lifecycle."""
import base64
import json
from pathlib import Path

import numpy as np
import pytest

from ipi.engine.forcefields import FFEval, FFSocket
from ipi.engine.outputs import CheckpointOutput
from ipi.engine.simulation import Simulation
from ipi.utils.io import NumpyEncoder


class RingField(FFSocket):
    """Synchronous socket stand-in; also serializes as an ordinary ffsocket."""
    poll = FFEval.poll
    _eval_thread = FFEval._eval_thread

    def __init__(self):
        super().__init__(name="ace", dopbc=False, threaded=False)
        self.calls = 0
        self.bad = None

    def start(self):
        pass

    def stop(self):
        pass

    def evaluate(self, request):
        self.calls += 1
        config = json.loads(base64.b64decode(self.pars["ring_polymer"]))
        q = request["pos"].reshape(config["nbeads"], -1)
        C = np.array(config["transform"])
        alpha = np.array(config["alpha"])
        d = q[:, :3] - q[:, :3].mean(axis=0)
        c = 1 + .2 * q[:, 3:]
        F = c * d
        w = C.T @ (alpha[:, None] * (C @ F))
        g = np.zeros_like(q)
        g[:, :3] = c * w - (c * w).mean(axis=0)
        g[:, 3:] = .2 * d * w
        jac = np.array([np.diag(x) for x in c])
        extras = dict(friction_ring_schema="ipi_friction_ring_v1",
                      friction_coupling_jacobian=jac.tolist(),
                      friction_gamma=np.einsum("bci,bcj->bij", jac, jac).tolist(),
                      friction_meta=dict(schema="ipi_friction_v1", active_atoms=[0],
                          channel_labels=["x", "y", "z"],
                          jacobian_units="ipi_internal_mass_scaled",
                          gamma_units="ipi_internal_mass_scaled"))
        energy = .5 * np.sum(alpha[:, None] * (C @ F)**2)
        force = -g.ravel()
        if self.bad == "schema":
            extras.pop("friction_ring_schema")
        elif self.bad == "force":
            force = force[:-1]
        elif self.bad == "nan":
            energy = np.nan
        request["result"] = [energy, force, np.zeros((3, 3)), extras]
        request["status"] = "Done"
        request._event_done.set()


def load(xml):
    def hook(sim):
        from types import MethodType
        field = FFSocket(name="ace", dopbc=False, threaded=False)
        field.calls, field.bad = 0, None
        for name in ("poll", "_eval_thread", "start", "stop", "evaluate"):
            setattr(field, name, MethodType(getattr(RingField, name), field))
        sim.fflist["ace"] = field
    return Simulation.load_from_xml(xml, post_init_hook=hook)


def create(folder, dt=.02, fixed="", spread=False):
    folder.mkdir(parents=True, exist_ok=True)
    driver = folder / "potential.py"
    driver.write_text('''import numpy as np
from ipi.pes.dummy import Dummy_driver
__DRIVER_CLASS__="Model"
class Model(Dummy_driver):
 def compute_structure(self,cell,pos):
  return .5*np.sum(pos**2),-pos,np.zeros((3,3)),""
''')
    q = np.array([[.1,.2,.3,1,1,1],[.2,.3,.4,1,1,1],[.3,.4,.5,1,1,1]])
    if spread:
        q[0,3] += .1
    p = np.tile([.2,.3,.1,0,0,0],(3,1))
    ap = np.tile([[0.,.2],[-.2,.3]],(3,1,1))
    def arr(tag, a):
        return f'<{tag} shape="{a.shape}">['+','.join(map(str,a.ravel()))+f']</{tag}>'
    fixed_xml = f'<fixatoms_dof>[{fixed}]</fixatoms_dof>' if fixed else ''
    xml = f'''<simulation verbosity="quiet" floatformat="%.17e">
    <prng><seed>421</seed></prng><total_steps>1000</total_steps>
    <ffsocket name="ace" mode="unix" pbc="false"><address>unused_test</address></ffsocket>
    <ffdirect name="physical" pbc="false"><pes>custom</pes><pes_path>{driver}</pes_path></ffdirect><output></output>
    <system><beads natoms="2" nbeads="3">{arr('q',q)}{arr('p',p)}
    <m shape="(2)">[1,2]</m><names shape="(2)">[H,Cu]</names></beads>
    <cell shape="(3,3)">[10,1,2,0,11,3,0,0,12]</cell>
    <ensemble><temperature>.01</temperature></ensemble>
    <forces><force forcefield="physical"/></forces>
    <motion mode="dynamics"><fixcom>False</fixcom>{fixed_xml}
    <dynamics mode="nve-f" splitting="symmetric-gle"><timestep>{dt}</timestep>
    <friction><forcefield>ace</forcefield><bath_mode>non-markovian</bath_mode>
    <variable_friction>True</variable_friction>{arr('Ap',ap)}
    <Lambda shape="(4,2)">[0,.5,1,.3,2,.1,10,0]</Lambda>
    </friction></dynamics></motion></system></simulation>'''
    return load(xml)


def advance(s):
    s.syslist[0].motion.step(step=s.step)
    s.step += 1


def state(s):
    sy = s.syslist[0]
    sy.motion.friction.bath._get_non_markovian_jacobian()
    return [sy.beads.q.copy(),sy.beads.p.copy(),sy.motion.friction.bath.s.copy(),
            np.array(sy.motion.friction.ediss)]


def test_force_cache_and_restart(tmp_path):
    s = create(tmp_path)
    fr = s.syslist[0].motion.friction
    b = s.syslist[0].beads
    provider = s.fflist["ace"]
    f = fr.force_mf.copy()
    q = b.q.copy()
    calls = provider.calls
    _ = fr.energy_mf, fr.coupling_jacobian, fr.gamma
    assert provider.calls == calls  # one evaluation supplies all four results
    for index in np.ndindex(q.shape):
        delta = np.zeros_like(q); delta[index] = 1e-5
        b.q = q + delta; ep = fr.energy_mf
        b.q = q - delta; em = fr.energy_mf
        np.testing.assert_allclose(f[index], -(ep-em)/2e-5, atol=1e-10)
    b.q = q
    for _ in range(7):
        advance(s)
    cp = CheckpointOutput(filename=str(tmp_path/'restart'),stride=1)
    cp.bind(s); cp.write(store=False)
    expected = []
    for _ in range(10):
        advance(s); expected.append(state(s))
    rng = json.dumps(s.prng.state,cls=NumpyEncoder)
    resumed = load((tmp_path/'restart').read_text())
    for ref in expected:
        advance(resumed)
        for actual,wanted in zip(state(resumed),ref):
            np.testing.assert_allclose(actual,wanted,atol=1e-13,rtol=1e-12)
    assert json.dumps(resumed.prng.state,cls=NumpyEncoder) == rng


@pytest.mark.parametrize("bad", ["schema", "force", "nan"])
def test_reject_bad_reply(tmp_path, bad):
    s = create(tmp_path)
    s.fflist["ace"].bad = bad
    with pytest.raises(ValueError):
        _ = s.syslist[0].motion.friction.energy_mf
