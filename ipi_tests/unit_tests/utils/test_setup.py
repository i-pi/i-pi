"""Deals with testing the automated driver building system.

Note that this will only run if you have Python version 3 or later.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import os
from ipi import install_driver


def test_driver_base():
    """Tests that running the driver install works"""

    install_driver()


def test_driver_noroot():
    """Tests what happens if IPI_ROOT is not set"""

    if "IPI_ROOT" in os.environ:
        ipi_root = os.environ["IPI_ROOT"]
        del os.environ["IPI_ROOT"]
        install_driver()
        os.environ["IPI_ROOT"] = ipi_root


def test_driver_forcebuild(tmp_path):
    """Tests that the driver is compiled and installed in the requested place"""

    driver = tmp_path / "i-pi-driver"
    install_driver(force_install=True, path=driver)
    assert os.access(driver, os.X_OK)
