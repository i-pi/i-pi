import subprocess as sp
from pathlib import Path
import pytest
import numpy as np
import time
import os
import shutil
from tempfile import TemporaryDirectory
import tempfile


def get_info_test(parent):
    folders = [x[0] for x in os.walk(parent)]
    reg_tests = list()

    for ff in folders:
        if os.path.isfile(Path(ff) / "input.xml"):
             driver_info = list()
             try:
                with open(Path(ff) / "driver.txt") as f:
                    ncount =0
                    while True:
                          line = f.readline()
                          if not line:
                              break
                          else:
                              driver, address, port , mode = line.split()             
                              driver_info.append((driver, address, port, mode))
                              ncount +=1
                    reg_tests.append((ff, driver_info))
                    if ncount == 0:
                         raise ValueError('driver.txt is empty')
             except FileNotFoundError:
                  raise FileNotFoundError('A driver.txt file is needed for each reg_text to be executed')
             except ValueError:
                  raise ValueError('Please specify "model address      port mode" inside driver.txt')

    return reg_tests


class Runner(object):
    def __init__(self, parent):
        self.parent = parent

    def _run(self, cmd1, cmd2, cwd):
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
               
            self._check_error(ipi)
            self._check_properties(cwd)

        except sp.TimeoutExpired:
            raise RuntimeError(
                "Time is out. Aborted during {} test. \
              Error {}".format(
                    str(cwd), ipi.communicate(timeout=2)[0]
                )
            )

        except AssertionError:
            raise AssertionError("{}".format(str(cwd)))

        except FileNotFoundError:
            raise FileNotFoundError("{}".format(str(cwd)))

        except ValueError:
            raise ValueError("{}".format(str(cwd)))

    def _check_error(self, ipi):
        ipi_error = ipi.communicate(timeout=30)[1].decode("ascii")
        print(ipi_error)
        assert "" == ipi_error

    def _check_properties(self,cwd):
        try:
           ref_output = np.loadtxt(Path(cwd)/'ref_simulation.out')
        except IOError:
           raise ('Please provide a refence properties output named "ref_simulation.out"')
        except ValueError:
            raise ValueError("Please check ref_simulation.out in {}".format(str(self.parent)))
   
        test_output = np.loadtxt(self.tmp_dir/'simulation.out')
 
        np.testing.assert_array_equal( test_output, ref_output)

