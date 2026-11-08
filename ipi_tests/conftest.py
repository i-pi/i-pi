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
    parser.addoption(
        "--ipi-mpi",
        action="store_true",
        help="run the end-to-end tests of the MPI transport, which need mpi4py, "
        "an MPI runtime and the driver built with `make MPI=1`",
    )
    parser.addoption(
        "--ipi-metatomic-model",
        default="",
        help="absolute path of the PET-MAD-xs v1.6.0 model of the metatomic tests "
        "(by default, it is exported with upet)",
    )
    parser.addoption(
        "--ipi-petjax-checkpoint",
        default="",
        help="absolute path of the pet-jax checkpoint of PET-MAD-xs v1.6.0, which "
        "the petjax tests compare with the metatomic model",
    )
    parser.addoption(
        "--ipi-fetch-driver",
        action="store_true",
        help="test fetching the Fortran driver from GitHub and compiling it",
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
    test_tools.ipi_mpi = config.getoption("--ipi-mpi")

    if test_tools.ipi_backend in ("torch", "jax"):
        # the reference outputs are generated with numpy: the Gaussian arrays
        # are drawn from the numpy generator also with torch and JAX, so that
        # runs in the same process follow the reference trajectories
        from ipi.utils.prng import Random

        for name, method in test_tools.numpy_draws().items():
            if name.endswith(test_tools.ipi_backend):
                _device_random[name] = getattr(Random, name)
                setattr(Random, name, method)

    if test_tools.ipi_backend == "torch":
        import torch

        # torch is slower than numpy for the small systems of the tests
        test_tools.TIMEOUT = 4 * test_tools.TIMEOUT

        # implicit conversions of tensors to numpy fail on GPUs: they are made to
        # fail also on the CPU, so that tests on the CPU detect them
        torch.Tensor.__array__ = _no_implicit_numpy

    if test_tools.ipi_backend == "jax":
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
