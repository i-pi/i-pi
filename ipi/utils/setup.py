""" Utility functions to install command-line tools.

These utilities facilitate the installation of the Fortran driver 
when it is not already present.
"""

import os
import shutil
import tempfile
import subprocess
from ipi.utils.messages import info, verbosity

__all__ = ['install_driver']

def get_ipi_path():
    """ Guesses the path where i-PI is installed """
    path = os.getenv('IPI_ROOT')
    if path is None:
        path = os.path.normpath(os.path.join(
            os.path.abspath(__file__),
            '../../../'
        ))

    return path

def install_driver():
    print(get_ipi_path())

    temp_dir = tempfile.mkdtemp()

    info(f"Temporary directory created: {temp_dir}", verbosity.medium)

    # fetches the drivers/f90 folder from the i-pi repo
    try:
        subprocess.run(['git', 'init'], cwd=temp_dir)
        subprocess.run(['git', 'remote', 'add', 'origin', 'https://github.com/i-pi/i-pi.git'], cwd=temp_dir)
        subprocess.run(['git', 'config', 'core.sparseCheckout', 'true'], cwd=temp_dir)
        with open(os.path.join(temp_dir, '.git', 'info', 'sparse-checkout'), 'w') as f:
            f.write('drivers/f90\n')
        subprocess.run(['git', 'pull', '--depth=1', 'origin', 'main'], cwd=temp_dir)
    except:
        info("Failed to fetch the drivers folder from i-PI github repository")
        raise

    # compiles the driver
    try:
        subprocess.run(['make', '-j', '1', 'driver.x'], cwd=os.path.join(temp_dir, 'drivers/f90/'))
    except:
        info("Failed to compile the FORTRAN driver. Make sure you have `gfortran` and `make` available")
        raise

    # moves the file to the bin/ folder
    try:
        ipi_bin = os.path.join(get_ipi_path(), 'bin')

        shutil.copy(os.path.join(temp_dir, 'drivers/f90/driver.x'),
                    os.path.join(ipi_bin, 'i-pi-driver')
        )
    except:
        info("Failed to copy the driver to the bin folder")
        raise


