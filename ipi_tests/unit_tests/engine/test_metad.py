"""Metadynamics must invalidate all outputs of the updated bias component."""

from types import SimpleNamespace
import numpy as np
import pytest
from ipi.engine.smotion.metad import MetaDyn
from ipi.utils.depend import depend_value


@pytest.mark.parametrize("work", [0.0, -0.25])
@pytest.mark.parametrize("weight", [0.0, 1.0])
def test_update_refreshes_only_active_bias(work, weight):
    """Zero work can still change force; refresh lazily and preserve other caches."""
    state = {"updated": False, "evaluations": 0, "updates": 0}

    def evaluate():
        state["evaluations"] += 1
        changed = state["updated"]
        return [
            1.0 - work if changed else 1.0,
            np.full(3, 2.0 if changed else 0.0),
            np.eye(3) * (3.0 if changed else 0.0),
            {"step": int(changed)},
        ]

    def update(pos, cell):
        state["updates"] += 1
        state["updated"] = True
        return work

    active = depend_value(name="active", func=evaluate)
    other = depend_value(name="other", func=lambda: [4.0])
    active.get()
    other.get()
    force = depend_value(
        name="force", func=lambda: active.get()[1], dependencies=[active]
    )
    force.get()
    component = lambda name, cache: SimpleNamespace(
        ffield=name, _forces=[SimpleNamespace(_ufvx=cache)]
    )
    system = SimpleNamespace(
        beads=SimpleNamespace(qc=np.zeros(3), nbeads=4),
        cell=SimpleNamespace(h=np.eye(3)),
        ensemble=SimpleNamespace(
            bcomp=[SimpleNamespace(ffield="bias")],
            bweights=[weight],
            eens=0.0,
            bias=SimpleNamespace(
                ff={"bias": SimpleNamespace(mtd_update=update)},
                mforces=[component("bias", active), component("other", other)],
            ),
        ),
    )
    motion = MetaDyn(metaff=["bias"])
    motion.syslist = [system]
    motion.step()
    assert state["evaluations"] == 1  # Invalidation does not evaluate eagerly.
    assert not other.tainted()
    if weight == 0:
        assert state["updates"] == 0
        assert not active.tainted()
        np.testing.assert_array_equal(force.get(), np.zeros(3))
    else:
        assert state["updates"] == 1
        assert active.tainted()
        np.testing.assert_array_equal(force.get(), np.full(3, 2.0))
        result = active.get()
        assert result[0] == 1.0 - work
        np.testing.assert_array_equal(result[2], np.eye(3) * 3.0)
        assert result[3] == {"step": 1}
        assert state["evaluations"] == 2  # Repeated reads reuse the fresh cache.
    assert system.ensemble.eens == (4 * work if weight else 0.0)
