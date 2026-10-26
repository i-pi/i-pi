"""Lattice and atom-identity invariants for the alloy KMC model."""

import numpy as np
from ipi.engine.motion.al6xxx_kmc import AlKMC


def test_fcc_neighbors_and_vacancy_exchange_preserve_species():
    kmc = AlKMC(
        mode="kmc",
        optimizer={"mode": "sd"},
        nstep=1,
        a0=4.0,
        ncell=4,
        nvac=1,
        nsi=2,
        nmg=2,
        neval=1,
        diffusion_barrier_al=0.1,
        diffusion_prefactor_al=2.0,
        diffusion_barrier_mg=-1,
        diffusion_prefactor_mg=-1,
        diffusion_barrier_si=-1,
        diffusion_prefactor_si=-1,
    )
    assert len(kmc.sites) == 64
    for i, neighbors in enumerate(kmc.neigh):
        assert len(set(neighbors)) == 12 and i not in neighbors
        for j in neighbors:
            assert i in kmc.neigh[j]
            displacement = kmc.dcell.array_pbc(kmc.sites[i] - kmc.sites[j])
            assert np.isclose(np.dot(displacement, displacement), 8.0)
    labels = np.array(list("SSMM" + "A" * 59 + "V"))
    indices = kmc.unique_idx(labels)
    np.testing.assert_array_equal(labels[indices], labels)
    changed = labels.copy()
    changed[4], changed[63] = changed[63], changed[4]
    mapping = kmc.unique_idx_match(labels, changed)
    moved = kmc.unique_idx(changed)[mapping[: kmc.natoms]]
    assert moved[4] == 63
    np.testing.assert_array_equal(
        np.delete(moved, 4), np.delete(indices[: kmc.natoms], 4)
    )
    assert kmc.barriers == {"A": 0.1, "M": 0.1, "S": 0.1}
