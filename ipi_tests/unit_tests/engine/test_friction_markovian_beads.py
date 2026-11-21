"""Reject beadwise Markovian dynamics outside the classical limit."""

from types import SimpleNamespace

import numpy as np
import pytest

from ipi.engine.friction import Friction, FrictionGLE


@pytest.mark.parametrize("variable", [False, True])
@pytest.mark.parametrize("mf", ["on", "off"])
@pytest.mark.parametrize("entry", ["bind", "step", "bath_bind"])
def test_multibead_markovian_is_rejected(variable, mf, entry):
    friction = Friction(
        variable_friction=variable, bath_mode="markovian", debug_mf_mode=mf,
    )
    beads = SimpleNamespace(nbeads=8, p=np.ones((8, 3)))
    friction.beads = beads
    motion = SimpleNamespace(
        beads=beads, nm=None, ensemble=SimpleNamespace(forces=None), prng=None,
    )
    with pytest.raises(ValueError, match="Markovian electronic friction requires nbeads=1"):
        if entry == "bind":
            friction.bind(motion)
        elif entry == "step":
            friction.step(0.1)
        else:
            FrictionGLE().bind(friction, None)
    np.testing.assert_array_equal(beads.p, np.ones((8, 3)))


@pytest.mark.parametrize("mode,nbeads", [("markovian", 1), ("non-markovian", 8), ("none", 8)])
def test_allowed_bath_bead_counts(mode, nbeads):
    friction = Friction(bath_mode=mode)
    friction.beads = SimpleNamespace(nbeads=nbeads)
    friction._validate_bath_bead_count()
