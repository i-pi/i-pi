#!/usr/bin/env python3

import pytest
from ....common import xyz_generator as xyz_gen
import ipi.utils.mathtools as mt
import ipi.utils.io.io_units as testing

from ipi.utils.units import Elements
import numpy as np
import numpy.testing as npt


test_data = [
    (
        1,
        1,
        "asdasd positions{angstrom}  100 aaa cell{angstrom} asdasd ",
        1.8897261,
        1.8897261,
    ),
    (
        5,
        1,
        "asdasd positions{angstrom}  100 aaa cell{angstrom} asdasd ",
        1.8897261,
        1.8897261,
    ),
    (
        1,
        1,
        "asdasd positions{atomic_units}  100 aaa cell{angstrom} asdasd ",
        1.0,
        1.8897261,
    ),
    (
        1,
        1,
        "asdasd positions{atomic_units}  100 aaa cell{atomic_units} asdasd ",
        1.0,
        1.0,
    ),
    (
        1,
        1,
        "asdasd positions{atomic_units}  100 aaa cell{meter} asdasd ",
        1.0,
        1.8897261e10,
    ),
]


@pytest.fixture(params=test_data)
def units_preparation(request):
    return request.param


@pytest.mark.skip(reason="This needs to be updated to match current code.")
def test_process_units_noobj(units_preparation):
    natoms, frames, comment, conver_xyz, conver_cell = units_preparation
    output = "noobj"

    filedesc, xyz, atom_names = xyz_gen.xyz_traj_filedesc(natoms, frames, comment)

    cell = mt.abc2h(1.0, 1.0, 1.0, np.pi / 2.0, np.pi / 2.0, np.pi / 2.0)

    masses = []
    _count = 0
    print("atom_names", atom_names)
    for _at in atom_names:
        print(Elements.mass(_at), _count, _at)
        masses.append(Elements.mass(_at))
        _count += 1

    masses = np.array(masses)
    type(masses)
    res = testing.process_units(
        comment,
        cell.copy(),
        xyz.copy(),
        np.array(atom_names).copy(),
        np.array(masses).copy(),
        output=output,
    )

    print(xyz, res["data"])
    npt.assert_array_almost_equal(res["data"], xyz * conver_xyz, 5)
    npt.assert_array_almost_equal(res["masses"], masses, 5)
    npt.assert_array_almost_equal(res["cell"], cell * conver_cell, 5)
    npt.assert_array_equal(res["names"], atom_names, 5)
    assert res["natoms"] == natoms


@pytest.mark.skip(reason="This needs to be updated to match current code.")
def test_process_units_object(units_preparation):
    natoms, frames, comment, conver_xyz, conver_cell = units_preparation
    output = "objects"

    filedesc, xyz, atom_names = xyz_gen.xyz_traj_filedesc(natoms, frames, comment)

    cell = mt.abc2h(1.0, 1.0, 1.0, np.pi / 2.0, np.pi / 2.0, np.pi / 2.0)

    masses = []
    _count = 0

    for _at in atom_names:
        print(Elements.mass(_at), _count, _at)
        masses.append(Elements.mass(_at))
        _count += 1

    masses = np.array(masses)

    res = testing.process_units(
        comment,
        cell.copy(),
        xyz.copy(),
        np.array(atom_names).copy(),
        np.array(masses).copy(),
        output=output,
    )

    npt.assert_array_almost_equal(res["atoms"].q, xyz * conver_xyz, 5)
    npt.assert_array_almost_equal(res["atoms"].m, masses, 5)
    npt.assert_array_almost_equal(res["cell"].h, cell * conver_cell, 5)
    npt.assert_array_equal(res["atoms"].names, atom_names, 5)
