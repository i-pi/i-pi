import subprocess as sp
from pathlib import Path
import pytest
import time
import os
import shutil
import tempfile


class Runner(object):
    def __init__(self, parent):
        self.parent = parent

    def _run(self, cmd1, cmd2, cwd):

        try:
            tmp_dir = Path(tempfile.mkdtemp())
            shutil.copytree(self.parent / cwd, tmp_dir / cwd)

            ipi = sp.Popen(
                cmd1, cwd=(tmp_dir / cwd), shell=True, stdout=sp.PIPE, stderr=sp.PIPE
            )
            time.sleep(3)
            driver = sp.Popen(cmd2, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
            self._check_error(ipi, driver)
            shutil.rmtree(tmp_dir)

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

    def _check_error(self, ipi, driver):
        assert "" == ipi.communicate(timeout=30)[1].decode("ascii")
