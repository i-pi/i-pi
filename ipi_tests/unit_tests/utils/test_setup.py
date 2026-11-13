"""Tests the installation of the Fortran driver (`ipi.install_driver`)."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2026 i-PI developers
# See the "licenses" directory for full license information.

import os
from pathlib import Path

import pytest

import ipi
from ipi import install_driver
from ipi.utils.setup import get_ipi_path

# the source tree of the tests, which holds drivers/f90
ROOT = Path(__file__).resolve().parents[3]


def test_existing_driver_is_kept(tmp_path, monkeypatch):
    """A driver that is already on PATH is not compiled again"""

    driver = tmp_path / "i-pi-driver"
    driver.write_text("#!/bin/sh\n")
    driver.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
    install_driver()
    assert driver.read_text() == "#!/bin/sh\n"


def test_root_is_guessed_from_the_package(monkeypatch):
    """Without IPI_ROOT, the drivers are looked for next to the ipi package"""

    monkeypatch.delenv("IPI_ROOT", raising=False)
    assert Path(get_ipi_path()) == Path(ipi.__file__).resolve().parents[1]


def test_driver_is_compiled_from_the_tree(tmp_path, monkeypatch):
    """The driver is compiled from the drivers/f90 folder of IPI_ROOT, and
    installed in the requested place"""

    monkeypatch.setenv("IPI_ROOT", str(ROOT))
    driver = tmp_path / "i-pi-driver"
    install_driver(force_install=True, path=driver)
    assert os.access(driver, os.X_OK)


def test_driver_is_fetched_from_github(tmp_path, monkeypatch, request):
    """Without a drivers/f90 folder, the driver is fetched from GitHub"""

    if not request.config.getoption("--ipi-fetch-driver"):
        pytest.skip("run with --ipi-fetch-driver to fetch the driver from GitHub")
    monkeypatch.setenv("IPI_ROOT", str(tmp_path / "empty"))
    driver = tmp_path / "i-pi-driver"
    install_driver(force_install=True, path=driver)
    assert os.access(driver, os.X_OK)
