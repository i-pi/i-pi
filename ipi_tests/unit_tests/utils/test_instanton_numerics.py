"""Spring matrices, path doubling and interpolation in instanton utilities."""

from types import SimpleNamespace
import numpy as np
import pytest
from scipy.linalg import block_diag
from ipi.engine.beads import Beads
from ipi.utils.instools import (
    banded_hessian,
    invmul_banded,
    red2comp,
    get_imvector,
    Fix,
)
from ipi.utils.tools.instanton_postproc import get_double, get_rp_freq, spring_pot
from ipi.utils.tools.instanton_interpolation import interpolate_instanton
from ipi.utils.io import read_file
from ipi.utils.units import unit_to_internal


@pytest.mark.parametrize("nbeads", [1, 3, 5])
@pytest.mark.parametrize("mass_weighted", [False, True])
def test_banded_newton_solution_matches_dense_spring_matrix(nbeads, mass_weighted):
    beads = Beads(1, nbeads)
    beads.m[:] = 2.0
    block = np.array([[3.0, 0.2, 0.1], [0.2, 4.0, 0.3], [0.1, 0.3, 5.0]])
    reduced = np.tile(block, (1, nbeads))
    mapper = SimpleNamespace(
        fix=SimpleNamespace(fixbeads=beads), coef=np.ones(nbeads + 1), omegan=0.7
    )
    band = banded_hessian(reduced, mapper, masses=mass_weighted, shift=0.01)
    dense = block_diag(*[block] * nbeads) + 0.01 * np.eye(3 * nbeads)
    for b in range(nbeads - 1):
        for d in range(3):
            bond = np.zeros(3 * nbeads)
            bond[3 * b + d] = 1
            bond[3 * (b + 1) + d] = -1
            dense += 0.7**2 * (2 if mass_weighted else 1) * np.outer(bond, bond)
    rhs = np.arange(3 * nbeads) + 1.0
    for posdef in [False, True]:
        answer = invmul_banded(band, rhs, posdef)
        np.testing.assert_allclose(dense @ answer, rhs, atol=1e-12)
    np.testing.assert_allclose(
        red2comp(reduced, nbeads, 1), block_diag(*[block] * nbeads)
    )


@pytest.mark.parametrize("nbeads", [1, 4, 7])
@pytest.mark.parametrize("mode", ["rate", "splitting"])
def test_ring_polymer_frequencies_match_dense_chain(nbeads, mode):
    omega2 = np.array([0.3, 0.8, 1.2])
    temp = 0.4
    if mode == "rate":
        lap = np.zeros((nbeads, nbeads))
        for b in range(nbeads):
            v = np.zeros(nbeads)
            v[b] += 1
            v[(b + 1) % nbeads] -= 1
            lap += np.outer(v, v)
    else:
        lap = 2 * np.eye(nbeads) - np.eye(nbeads, k=1) - np.eye(nbeads, k=-1)
    values = np.concatenate(
        [
            np.linalg.eigvalsh(lap * (nbeads * temp) ** 2 + w * np.eye(nbeads))
            for w in omega2
        ]
    )
    actual = get_rp_freq(omega2, nbeads, temp, mode)
    if mode == "rate":
        assert actual == pytest.approx(
            np.log(np.sqrt(values) / (nbeads * temp)).sum(), abs=1e-12
        )
    else:
        np.testing.assert_allclose(
            np.sort(actual), np.sqrt(np.sort(values)), atol=1e-12
        )


def test_double_path_preserves_reflection_and_spring_energy():
    q = np.arange(12.0).reshape(4, 3) / 10
    h = block_diag(*[np.eye(3) * (i + 1) for i in range(4)])
    full, n, hfull = get_double(q, 4, 1, h)
    assert n == 8
    np.testing.assert_array_equal(full[:4], full[4:][::-1])
    np.testing.assert_array_equal(
        hfull, block_diag(h, *[np.eye(3) * i for i in [4, 3, 2, 1]])
    )
    assert spring_pot(8, full, 0.2, np.ones_like(full) * 2) == pytest.approx(
        2 * spring_pot(4, q, 0.2, np.ones_like(q) * 2)
    )


def test_unstable_direction_and_invalid_stable_hessian():
    v = get_imvector(np.diag([-2.0, 3.0, 4.0]), np.array([2.0, 3.0, 4.0]))
    np.testing.assert_allclose(np.abs(v), [[1, 0, 0]])
    with pytest.raises(ValueError, match="positive"):
        get_imvector(np.eye(3), np.ones(3))
    with pytest.raises(ValueError, match="size"):
        get_imvector(np.eye(6), np.ones(3))


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="Interpolated Hessian output remains buffered after the function returns",
)
@pytest.mark.parametrize("new_count", [2, 4, 8])
def test_interpolation_of_constant_instanton(tmp_path, monkeypatch, new_count):
    monkeypatch.chdir(tmp_path)
    with open("path.xyz", "w") as f:
        for _ in range(4):
            f.write(
                "1\n# CELL(abcABC): 10 10 10 90 90 90 positions{angstrom} cell{angstrom}\nH .2 .3 .4\n"
            )
    block = np.array([[2.0, 0.1, 0.2], [0.1, 3.0, 0.3], [0.2, 0.3, 4.0]])
    np.savetxt("hessian", np.tile(block, (1, 4)).reshape(1, -1))
    interpolate_instanton("None", "path.xyz", "hessian", new_count, True)
    with open("new_instanton.xyz") as f:
        frames = [read_file("xyz", f)["atoms"].q.copy() for _ in range(new_count)]
        with pytest.raises(EOFError):
            read_file("xyz", f)
    np.testing.assert_allclose(
        frames,
        np.tile([0.2, 0.3, 0.4], (new_count, 1))
        * unit_to_internal("length", "angstrom", 1.0),
        atol=1e-10,
    )
    np.testing.assert_allclose(
        np.loadtxt("new_hessian.dat").reshape(3, -1),
        np.tile(block, (1, new_count)),
        atol=1e-12,
    )


@pytest.mark.parametrize(
    "kind",
    [
        1,
        2,
        pytest.param(
            3,
            marks=pytest.mark.xfail(
                strict=True,
                raises=AttributeError,
                reason="Fix.get_full_vector references nonexistent self.fix",
            ),
        ),
        pytest.param(
            4,
            marks=pytest.mark.xfail(
                strict=True,
                raises=IndexError,
                reason="Friction Hessian expansion does not reset inner indices",
            ),
        ),
        pytest.param(
            5,
            marks=pytest.mark.xfail(
                strict=True,
                raises=IndexError,
                reason="Friction matrix expansion does not reset inner index",
            ),
        ),
    ],
)
def test_fixed_atom_reduction_and_expansion_roundtrip(kind):
    beads = Beads(2, 3)
    beads.m[:] = [2, 3]
    fix = Fix(np.array([0, 1, 2]), beads)
    shapes = {1: (3, 6), 2: (6, 18), 3: (2, 18), 4: (3, 6, 6, 6), 5: (3, 6, 6)}
    a = np.arange(np.prod(shapes[kind]), dtype=float).reshape(shapes[kind])
    active = fix.get_active_vector(a, kind)
    restored = fix.get_full_vector(active, kind)
    np.testing.assert_array_equal(fix.get_active_vector(restored, kind), active)
    assert np.count_nonzero(restored) <= active.size
