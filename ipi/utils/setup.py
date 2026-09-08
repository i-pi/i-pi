"""Utility functions to install command-line tools.

These utilities facilitate the installation of the Fortran driver
when it is not already present.
"""

import os
import shutil
import tempfile
import subprocess
from ipi.utils.messages import info, warning, verbosity

__all__ = ["install_driver"]


def get_ipi_path():
    """Guesses the path where i-PI is installed"""

    path = os.getenv("IPI_ROOT")
    if path is None or path == "":
        path = os.path.normpath(os.path.join(os.path.abspath(__file__), "../../"))

    return path


def install_driver(force_install=False):
    """
    This utility function fetches the FORTRAN driver folder from the ipi repo and installs it.
    Requires a system with git, gfortran and make.
    """

    ipi_path = get_ipi_path()
    ipi_driver_path = shutil.which("i-pi-driver")
    if ipi_driver_path is None or force_install:
        import sysconfig

        # get potential installation locations for the driver
        script_dirs = [
            os.path.join(ipi_path, "bin"),  # preferred choice
            sysconfig.get_path("scripts"),  # system, conda or virtual environment
            sysconfig.get_path("scripts", sysconfig.get_preferred_scheme("user")),
        ]

        os_path = os.getenv("PATH")

        if os_path is None or os_path == "":
            path_dirs = []
        else:
            path_dirs = list(map(os.path.abspath, os_path.split(os.pathsep)))

        for directory in script_dirs:
            # check if the directory is on PATH and writable
            if directory in path_dirs and os.access(directory, os.W_OK):
                ipi_driver_path = os.path.join(directory, "i-pi-driver")
                break

        if ipi_driver_path is None:
            # install locally
            ipi_driver_path = os.path.join(os.getcwd(), "i-pi-driver")
    else:
        info(f"i-pi-driver is already present in {ipi_driver_path} ", verbosity.low)
        return

    build_dir = os.path.join(ipi_path, "drivers/f90")

    if not os.path.exists(build_dir):
        # fetches the drivers/f90 folder from the i-pi repo
        temp_dir = tempfile.mkdtemp()

        info(f"Temporary directory created: {temp_dir}", verbosity.medium)

        try:
            subprocess.run(["git", "init"], cwd=temp_dir)
            subprocess.run(
                [
                    "git",
                    "remote",
                    "add",
                    "origin",
                    "https://github.com/lab-cosmo/i-pi.git",
                ],
                cwd=temp_dir,
            )
            subprocess.run(
                ["git", "config", "core.sparseCheckout", "true"], cwd=temp_dir
            )
            with open(
                os.path.join(temp_dir, ".git", "info", "sparse-checkout"), "w"
            ) as f:
                f.write("drivers/f90\n")
            subprocess.run(["git", "pull", "--depth=1", "origin", "main"], cwd=temp_dir)
        except:
            warning(
                "Failed to fetch the drivers folder from i-PI github repository",
                verbosity.low,
            )
            raise

        build_dir = os.path.join(temp_dir, "drivers/f90")

    # compiles the driver
    try:
        subprocess.run(["make", "-j", "1", "driver.x"], cwd=build_dir)
    except:
        warning(
            "Failed to compile the FORTRAN driver. Make sure you have `gfortran` and `make` available",
            verbosity.low,
        )
        raise

    # moves the file to the bin/ folder
    try:
        shutil.copy(os.path.join(build_dir, "driver.x"), ipi_driver_path)
    except:
        warning("Failed to copy the driver to the bin folder", verbosity.low)
        raise
