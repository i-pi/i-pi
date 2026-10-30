"""Lattice and atom-identity invariants for the alloy KMC model."""

import numpy as np
import pytest

import ipi.utils.array_backend as ab
from ipi.engine.motion.al6xxx_kmc import AlKMC

# a 4x4x4 FCC lattice (a0 = 7.6 bohr) with one vacancy, and zero forces
KMC_INPUT = """<simulation verbosity="quiet" threading="False">
  <total_steps>10</total_steps><prng><seed>7</seed></prng>
  <ffdirect name="ff"><pes>dummy</pes></ffdirect>
  <system>
    <initialize nbeads="1"><file mode="xyz">init.xyz</file></initialize>
    <forces><force forcefield="ff"/></forces>
    <ensemble><temperature units="kelvin">300</temperature></ensemble>
    <motion mode="al-kmc">
      <al6xxx_kmc mode="rfkmc">
        <optimizer mode="sd"/>
        <nstep>1</nstep><a0>7.6</a0><ncell>4</ncell>
        <nvac>1</nvac><nsi>2</nsi><nmg>2</nmg><neval>2</neval>
        <diffusion_barrier_al>0.01</diffusion_barrier_al>
        <diffusion_prefactor_al>2.0</diffusion_prefactor_al>
      </al6xxx_kmc>
    </motion>
  </system>
</simulation>"""


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


@pytest.mark.parametrize("backend", ["numpy", "torch"])
def test_kmc_step_caches_host_values(simulation, backend, request):
    """A KMC step runs with the arrays of either backend, and the cached energies
    and structures, which are saved with pickle, are floats and numpy arrays"""
    if backend == "torch":
        request.getfixturevalue("strict_torch")
    # one particle per lattice site, the vacancy last
    names = ["Si"] * 2 + ["Mg"] * 2 + ["Al"] * 59 + ["V"]
    with open("init.xyz", "w") as f:
        f.write("64\n# CELL(abcABC): 30.4 30.4 30.4 90 90 90 ")
        f.write("cell{atomic_unit} positions{atomic_unit}\n")
        f.writelines(f"{name} 0 0 0\n" for name in names)
    ab.set_backend(backend)
    try:
        with simulation(KMC_INPUT) as sim:
            kmc = sim.syslist[0].motion
            sim.run(steps=1, write_outputs=False)
            assert len(kmc.ecache) > 1
            assert all(isinstance(e, float) for e in kmc.ecache.values())
            assert all(type(q) is np.ndarray for q in kmc.qcache.values())
    finally:
        ab.set_backend("numpy")
