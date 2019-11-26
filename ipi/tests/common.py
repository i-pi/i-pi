"""Common helper functions for running the tests."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import glob
import os
import subprocess
import shutil
import tempfile


def local(file=None):
    """Returns local folder of the tests directory.

    Args:
        - file: Append file to the local folder
    """
    if file is None:
        return os.path.dirname(__file__)
    else:
        return os.path.join(os.path.dirname(__file__), file)


class SimulationTest(object):

    """Simple class used to test various aspects of the simulation.

    Can be used to run an example given the location of an xml
    input file and the location of a suitable driver code.

    Attributes:
       finput: The name of the xml input file
       folder_input: A string giving the directory the input file is held in.
       fdriver: The location of a driver code.
       cwd: Current working directory.
       tmpdir: A temporary directory to run the simulation in.
    """

    def __init__(self, input, driver):
        """Initializes SimulationTest.

        Args:
           input: The name of the xml input file.
           driver: The location of the driver code.
        """

        self.finput = input
        self.folder_input = os.path.dirname(input)
        self.fdriver = driver
        self.cwd = os.getcwd()
        self.tmpdir = tempfile.mkdtemp()

        # Copy needed files to tmpdir
        for src in glob.glob("%s/*" % self.folder_input):
            shutil.copy(src, self.tmpdir)

        os.chdir(self.tmpdir)

    def __del__(self):
        """Cleans the temporary directory once the simulation is over."""

        os.chdir(self.cwd)
        shutil.rmtree(self.tmpdir)

    def run(self):
        """Runs the simulation."""

        # Run driver
        p = subprocess.Popen("echo running %s" % self.fdriver, shell=True)

        # Start simulation
        # TODO
        print(subprocess.check_output("ls", shell=True))
        print(subprocess.check_output("pwd", shell=True))

        # wait for driver to finish
        p.wait()
