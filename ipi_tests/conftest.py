"""pytest options and fixtures shared by the i-PI test suites."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import pytest

# methods of Random that draw on the device, replaced with --ipi-backend torch
# or jax
_device_random = {}


def pytest_addoption(parser):
    parser.addoption(
        "--ipi-backend",
        default="numpy",
        choices=["numpy", "torch", "jax"],
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

    if test_tools.ipi_backend in ("torch", "jax"):
        # the reference outputs are generated with numpy: the Gaussian arrays
        # are drawn from the numpy generator also with torch and JAX, so that
        # runs in the same process follow the reference trajectories
        from ipi.utils import array_backend
        from ipi.utils.prng import Random

        def gvec(self, shape):
            return array_backend.xp.asarray(self.gvec_serial(shape))

    if test_tools.ipi_backend == "torch":
        import torch

        def gfill(self, out):
            out[...] = array_backend.xp.asarray(self.gvec_serial(tuple(out.shape)))

        _device_random.update(
            gvec_torch=Random.gvec_torch, gfill_torch=Random.gfill_torch
        )
        Random.gvec_torch = gvec
        Random.gfill_torch = gfill

        # torch is slower than numpy for the small systems of the tests
        test_tools.TIMEOUT = 4 * test_tools.TIMEOUT

        # implicit conversions of tensors to numpy fail on GPUs: they are made to
        # fail also on the CPU, so that tests on the CPU detect them
        torch.Tensor.__array__ = _no_implicit_numpy

    if test_tools.ipi_backend == "jax":
        _device_random.update(gvec_jax=Random.gvec_jax)
        Random.gvec_jax = gvec
        # eager JAX is much slower than numpy for the small systems of the tests
        test_tools.TIMEOUT = 10 * test_tools.TIMEOUT


@pytest.fixture
def device_random(monkeypatch):
    """Draws the Gaussian arrays with the generator of torch or JAX, on the
    device, also when they are drawn from the numpy generator (--ipi-backend
    torch or jax)."""

    from ipi.utils.prng import Random

    for name, method in _device_random.items():
        monkeypatch.setattr(Random, name, method)


@pytest.fixture
def strict_torch(monkeypatch):
    """Makes implicit conversions of torch tensors to numpy arrays fail also on
    the CPU, as they do for tensors on a GPU."""

    torch = pytest.importorskip("torch")
    monkeypatch.setattr(torch.Tensor, "__array__", _no_implicit_numpy)
