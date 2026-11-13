"""Fixtures shared by the unit tests."""

import pytest

import ipi.utils.array_backend as ab


@pytest.fixture
def select_backend():
    """Selects an array backend, skipping the test if its library is missing,
    and restores the numpy backend after the test."""

    def select(name, device="cpu"):
        if name != "numpy":
            pytest.importorskip(name)
        ab.set_backend(name, device=device)

    yield select
    ab.set_backend("numpy")
