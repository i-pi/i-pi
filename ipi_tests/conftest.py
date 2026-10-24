"""pytest options and fixtures shared by the i-PI test suites."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--ipi-backend",
        default="numpy",
        choices=["numpy", "torch"],
        help="array backend used to run i-PI in the regression tests and examples",
    )
    parser.addoption(
        "--ipi-device",
        default="cpu",
        help="device used by the array backend (e.g. 'cuda')",
    )


def _no_implicit_numpy(self, *args, **kwargs):
    raise TypeError(
        "implicit conversion of a torch tensor to a numpy array: "
        "use ipi.utils.array_backend.to_numpy at the boundaries"
    )


def pytest_configure(config):
    from ipi_tests import test_tools

    test_tools.ipi_backend = config.getoption("--ipi-backend")
    test_tools.ipi_device = config.getoption("--ipi-device")

    if test_tools.ipi_backend == "torch":
        # the reference outputs are generated with numpy: the Gaussian arrays
        # are drawn from the numpy generator also with torch, so that runs
        # in the same process follow the reference trajectories
        import torch
        from ipi.utils import array_backend
        from ipi.utils.prng import Random

        def gvec(self, shape):
            return array_backend.xp.asarray(self.gvec_serial(shape))

        def gfill(self, out):
            out[...] = array_backend.xp.asarray(self.gvec_serial(tuple(out.shape)))

        Random.gvec_torch = gvec
        Random.gfill_torch = gfill

        # torch is slower than numpy for the small systems of the tests
        test_tools.TIMEOUT = 4 * test_tools.TIMEOUT

        # implicit conversions of tensors to numpy fail on GPUs: they are made to
        # fail also on the CPU, so that tests on the CPU detect them
        torch.Tensor.__array__ = _no_implicit_numpy


@pytest.fixture
def strict_torch(monkeypatch):
    """Makes implicit conversions of torch tensors to numpy arrays fail also on
    the CPU, as they do for tensors on a GPU."""

    torch = pytest.importorskip("torch")
    monkeypatch.setattr(torch.Tensor, "__array__", _no_implicit_numpy)
