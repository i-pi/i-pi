#!/usr/bin/env python3

"""See [documentation](https://janberges.github.io/elphmod/modules/md.html)."""

import elphmod

el = elphmod.el.Model("model", rydberg=True)
ph = elphmod.ph.Model("model.ifc", divide_mass=False)
elph = elphmod.elph.Model("model.epmatwp", "model.wigner", el, ph, divide_mass=False)

driver = elphmod.md.Driver(
    elph,
    nk=(12, 12),
    nq=(2, 2),
    supercell=(9, 9),
    kT=0.02,
    f=elphmod.occupations.marzari_vanderbilt,
    n=1.0,
)

driver.kT = 0.005
driver.f = elphmod.occupations.fermi_dirac
driver.random_displacements()

driver.save("driver.pickle")
driver.to_xyz("driver.xyz")
