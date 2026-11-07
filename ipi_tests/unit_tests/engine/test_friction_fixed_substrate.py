import json
import numpy as np
import pytest
from ipi.engine.outputs import CheckpointOutput
from ipi.utils.io import NumpyEncoder
from ipi_tests.unit_tests.engine.test_friction_external_meanfield import create, advance, state, load

def test_fixed_restart(tmp_path):
 s=create(tmp_path/'r',fixed='3,4,5');sy=s.syslist[0];fr=sy.motion.friction
 assert np.max(abs(fr.force_mf[:,3:]))>1e-8 # Nonzero reaction force must be excluded.
 for _ in range(5):advance(s)
 np.testing.assert_array_equal(sy.beads.q[:,3:],np.ones((3,3)));np.testing.assert_array_equal(sy.beads.p[:,3:],np.zeros((3,3)))
 cp=CheckpointOutput(filename=str(tmp_path/'restart'),stride=1);cp.bind(s);cp.write(store=False);expected=[]
 for _ in range(6):advance(s);expected.append(state(s))
 rng=json.dumps(s.prng.state,cls=NumpyEncoder);r=load((tmp_path/'restart').read_text())
 for ref in expected:
  advance(r)
  for a,b in zip(state(r),ref):np.testing.assert_allclose(a,b,rtol=1e-12,atol=1e-13)
 assert json.dumps(r.prng.state,cls=NumpyEncoder)==rng

@pytest.mark.parametrize('fixed,spread',[('0,1,2',False),('3,4,5',True)])
def test_reject_bad_constraints(tmp_path,fixed,spread):
 with pytest.raises(ValueError):
  s=create(tmp_path/'bad',fixed=fixed,spread=spread)
  advance(s)
def test_fixed_energy_order(tmp_path):
 errors=[]
 for dt in [.02,.01]:
  s=create(tmp_path/str(dt),dt,fixed='3,4,5');sy=s.syslist[0];fr=sy.motion.friction;fr.bath._get_non_markovian_jacobian();fr.bath.os_step=lambda _:None
  def energy():return float(sy.forces.pot+sy.nm.kin+sy.nm.vspring+fr.energy_mf+.5*np.sum(fr.bath.s**2))
  values=[energy()]
  for _ in range(round(1/dt)):advance(s);values.append(energy())
  errors.append(max(abs(np.array(values)-values[0])))
 assert 3.5<errors[0]/errors[1]<4.5,errors
