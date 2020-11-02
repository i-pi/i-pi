import subprocess as sp
import sys
from pathlib import Path
import pytest
import numpy as np
import time
import glob
import os
import shutil
from tempfile import TemporaryDirectory
import tempfile


def get_info_test(parent):
    """This function recursively searches for examples
    and check for the presence of the required additional info
    required (i.e. driver.txt)
    """
    folders = [x[0] for x in os.walk(parent)]
    reg_tests = list()

    for ff in folders:
        if os.path.isfile(Path(ff) / "input.xml"):
            driver_info = list()
            try:
                with open(Path(ff) / "driver.txt") as f:
                    ncount = 0
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        else:
                            driver, address, port, mode = line.split()
                            driver_info.append((driver, address, port, mode))
                            ncount += 1
                    reg_tests.append((ff, driver_info))
                    if ncount == 0:
                        raise ValueError("driver.txt is empty")
            except FileNotFoundError:
                raise FileNotFoundError(
                    "({}) An input.xml file was found but a driver.txt not".format(ff)
                )
            except ValueError:
                raise ValueError(
                    'Please specify "model address      port mode" inside driver.txt'
                )

    return reg_tests


class Runner(object):
    """This class handles the modification the creation of tmp directories,
    the i-pi call, the driver call, and finally
    it checks that the generated outputs are the expected ones.
    """

    def __init__(
        self, parent, check_errors=True, check_main_output=True, check_xyz_output=False
    ):
        """ Store parent directory and commands to call i-pi and driver """
        self.parent = parent
        self.check_error = check_errors
        self.check_main_output = check_main_output
        self.check_xyz_output = check_xyz_output

    def _run(self, cmd1, cmd2, cwd):
        """This function tries to run the example in a tmp folder and
        checks if ipi has ended without error.
        After that the output is checked against a reference
        arguments:
            cwd: folder where all the original examples files are stored
            cmd1: command to call i-pi
            cmd2: list of commands to call drivers
        """

        try:
            for filename in glob.glob("/tmp/ipi_*"):
                os.remove(filename)
        except:
            pass

        try:
            self.tmp_dir = Path(tempfile.mkdtemp())
            files = os.listdir(self.parent / cwd)
            for f in files:
                shutil.copy(self.parent / cwd / f, self.tmp_dir)

            ipi = sp.Popen(
                cmd1, cwd=(self.tmp_dir), shell=True, stdout=sp.PIPE, stderr=sp.PIPE
            )
            time.sleep(2)
            driver = list()
            for cmd in cmd2:
                driver.append(sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE))

            if self.check_error:
                self._check_error(ipi)
            if self.check_main_output:
                self._check_main_output(cwd)
            if self.check_xyz_output:
                self._check_xyz_output(cwd)

        except sp.TimeoutExpired:
            raise RuntimeError(
                "Time is out. Aborted during {} test. \
              Error {}".format(
                    str(cwd), ipi.communicate(timeout=2)[0]
                )
            )

        except FileNotFoundError:
            raise FileNotFoundError("{}".format(str(cwd)))

        except ValueError:
            raise ValueError("{}".format(str(cwd)))

    def _check_error(self, ipi):
        """ This function checks if ipi has exited with errors"""
        ipi_error = ipi.communicate(timeout=30)[1].decode("ascii")
        print(ipi_error)
        assert "" == ipi_error

    def _check_main_output(self, cwd):
        """ This function checks if the simulation.out is 'all_close' to the reference file provided"""

        try:
            ref_output = np.loadtxt(Path(cwd) / "ref_simulation.out")
        except IOError:
            raise IOError(
                'Please provide a reference properties output named "ref_simulation.out"'
            )
        except ValueError:
            raise ValueError(
                "Please check ref_simulation.out in {}".format(str(self.parent))
            )

        test_output = np.loadtxt(self.tmp_dir / "simulation.out")
        np.testing.assert_allclose(test_output, ref_output)

    def _check_xyz_output(self, cwd):

        raise NotImplementedError
