#!/usr/bin/env python3
"""Run or create reference for regtest i-pi.

Todo:
    * Add support for .gzip files.
    * Send kill() if terminate() is not enough!

Dependancy:
    * numpy
    * colorama (optional - colored output)

The structure of the tree that contains the test must be as follow:

    `-- inputfiles_path    --> Specified in the `config files`
    |-- test1              --> Will be active if present in the `input_files`
    |   |-- something.xml  --> **Only** a single xml file (ipi input)
    |   |-- input_driver   --> Input files for the driver (see below)
    |   `-- regtest-ref    --> Folder containing output files to compare with.
    |-- test_2
    |   |-- something.xml
    |   |-- input_driver
    |   `-- output


The input must contains a specific comment that specifies which
files are needed by the driver and which command should be used to run the
driver. The syntax of the lines is quite strict: only spaces can be different.

To specify which files the driver needs:

<!--REGTEST
DEPENDENCIES  h5o2.dms4B.coeff.com.dat h5o2.pes4B.coeff.dat h5o2+.xyz
COMMAND(10)    i-pi-driver -u -h REGTEST_SOCKET -m zundel
COMMAND       i-pi-driver -u -h REGTEST_SOCKET -m zundel
ENDREGTEST-->

`DEPENDECIES`: are the file needed to run the test.
`COMMAND(n)`: command used to run the driver. The `n` are the number of
    instance of the COMMAND to be created. In the example there will be a
    total of 11 instance of the driver running at once.

This script will also try to take care of avoiding socket overlap. The most
most importance think in this regard is that, in case the name of the socket
compare more than once in the "command line" above, only the first one will be
considered. For example: using the command "i-pi-driver -u -m zundel -h zundel"
would not work.


# How the script works

All the tasks to be performed on a single test are grouped into the class Test.
This class inherit from threading.Thread then everything which is into the run
method can be perfomed in parallel. Since it could be needed to perform
different tasks on the same Test obj, a method "copy" allows to create a new
thread keeping all the properties of the previous thread.

The parallelism is based on multithreading and subprocess. The architecture has
has been tought to be as flexible as possible.


            +------------ Main Program --------------+
            v                                        v
   +-----------------+                      +-----------------+
   |                 |                      |                 |
   |    QUEUE_ALL    |               +----->|    QUEUE_COM    |
   |                 |               |      |                 |
   +--------+--------+               |      +--------+--------+
            |                        |               |
            |                        |               |
            |                        |               |
            +-----> Run TEST 1 ------+               |
            |                        |               |
            +-----> Run TEST 2 ------+               |
            |                        |               +---> Run Compare/Copy
            +-----> Run TEST 3 ------+                           ^
            |                        |                           |
            +-----> Run TEST n ------+                           |
                       ^                                         |
                       |                                         |
                       +--------------------+--------------------+
                                            |
                                            |
                                         Threads


This approach creates two queues: the tasks queued on QUEUE_ALL will be run
in parallel while the tasks queued on QUEUE_COM will be run in serial.

"""


import argparse
import os
import queue
import re
import shutil
import shlex
import subprocess as sbps
import sys
import threading
import time
import xml.etree.ElementTree as etree

import numpy as np
import numpy.testing as npt

# pylint: disable=import-error
from ipi.engine.outputs import CheckpointOutput, PropertyOutput, TrajectoryOutput
from ipi.engine.properties import getkey
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml

# pylint: enable=import-error

RED = YELLOW = GREEN = RESET = ""
try:
    from colorama import Fore  # pylint: disable=wrong-import-position

    RED = Fore.RED
    YELLOW = Fore.YELLOW
    GREEN = Fore.GREEN
    RESET = Fore.RESET
    INFO = Fore.LIGHTRED_EX
except ImportError:
    RED = ""
    YELLOW = ""
    GREEN = ""
    RESET = ""
    INFO = ""


# Hardcoded settings ####
TIMEOUT_DRIVER = 600  # Maximum time the driver are allowded to run
TIMEOUT_IPI = 30  # Maximum time to wait after the driver are done
IPI_WAITING_TIME = 5  # Time to wait after i-pi has been started
#


# Compile them only once! pylint: disable=anomalous-backslash-in-string
REGTEST_STRING_RGX = re.compile(r"<!--\s*REGTEST\s+([\s\w\W]*)\s+ENDREGTEST\s*-->")
# pylint: enable=anomalous-backslash-in-string

try:
    IPI_ROOT_FOLDER = os.environ["IPI_ROOT"]
except KeyError:
    IPI_ROOT_FOLDER = os.path.abspath("../")

QUEUE_ALL = queue.Queue()  # Queue to access the "run"
QUEUE_COM = queue.Queue()  # Queue to access the "compare"


def main():
    """Manage the main process."""

    root_test_folder = _parser()["root_test_folder"]
    tests_list = _parser()["tests"]

    # Create the root folder of the run
    root_run = _parser()["root_run_folder"]
    create_dir(root_run)

    # If no --add-test has been specified, search for tests in all directories
    # +within the root_test_folder.
    if len(tests_list) == 0:
        tests_list = [""]

    # Retrieve the test list and build the QUEUE_ALL
    test_list = _build_test_index(root_test_folder, tests_list)
    index = _parser()["index"]
    for test in test_list:
        test_obj = Test(index=index, name=test[0], path=test[1], root_run=root_run)
        QUEUE_ALL.put(test_obj)
        index += 10

    print("Starting tests")
    running_test = []
    running_com = []
    if int(_parser()["nproc"]) > 1:
        if (
            answer_is_y(
                '"!W! REGTEST is not thread-safe. Some tests could crash. Do you want to continue (y/n)?'
            )
            is False
        ):
            sys.exit()
    try:
        while True:
            if len(running_test) < _parser()["nproc"]:
                try:
                    running_test.append(QUEUE_ALL.get_nowait())
                except queue.Empty:
                    pass
                else:
                    running_test[-1].generate_output = True
                    running_test[-1].start()

            for _ii, thr in enumerate(running_test):
                thr.join(0.5)
                if not thr.is_alive():
                    thr = thr.copy()
                    thr.generate_output = False
                    QUEUE_COM.put(thr)
                    QUEUE_ALL.task_done()
                    del running_test[_ii]

            if len(running_com) < 1:
                try:
                    running_com.append(QUEUE_COM.get_nowait())
                except queue.Empty:
                    pass
                else:
                    running_com[-1].copy_reference = _parser()["create_reference"]
                    running_com[-1].compare_output = not _parser()["create_reference"]
                    running_com[-1].start()

            else:
                running_com[-1].join(0.5)
                if not running_com[-1].is_alive():
                    QUEUE_COM.task_done()
                    del running_com[-1]

            # print running_test, QUEUE_ALL.qsize(), QUEUE_COM.qsize(), running_com
            if (
                len(running_test) == 0
                and QUEUE_ALL.empty()
                and QUEUE_COM.empty()
                and len(running_com) == 0
            ):
                break

    except KeyboardInterrupt:
        for _thr in running_test:
            _thr.die = True

    print()


def _parser():
    """Parse the argument lists given as input.

    Return:
        A dictionary containing all the input options and argument.
    """
    parser = argparse.ArgumentParser(description="Regtests for i-PI.")
    parser.add_argument(
        "--maxtime",
        action="store",
        type=int,
        help=(
            "Wall time for the driver: this is useful when"
            "the driver is stuck. The default value is chosen"
            "based on the provided regtests."
        ),
        default=TIMEOUT_DRIVER,
        dest="maxtime",
    )
    parser.add_argument(
        "--add-test",
        action="append",
        type=str,
        help=(
            "Mark a test to be ran. A single test can be "
            'specified with "group:test". This option can be '
            "used as many times as you want."
        ),
        default=[],
        dest="tests",
    )
    parser.add_argument(
        "--tests-folder",
        action="store",
        type=str,
        default=os.path.join(IPI_ROOT_FOLDER, "examples"),
        help=("The folder where to search for tests."),
        dest="root_test_folder",
    )
    parser.add_argument(
        "--folder-run",
        action="store",
        type=str,
        default="regtest-run",
        help=(
            "Parent folder where the test will be run. "
            "If already existing all the content will "
            "be lost!"
        ),
        dest="root_run_folder",
    )
    parser.add_argument(
        "--starting-address",
        action="store",
        type=int,
        default=10000,
        help=("Number of the first port in case of inet " "socket."),
        dest="index",
    )
    parser.add_argument(
        "-np",
        "--nproc",
        action="store",
        type=int,
        default=1,
        help=("Number of concurrent test run at once."),
        dest="nproc",
    )
    parser.add_argument(
        "--create-reference",
        action="store_true",
        default=False,
        help=(
            "Only run the test, do not compare and at the "
            "end ask if you want to copy the output files "
            "in the reference test folder."
        ),
        dest="create_reference",
    )
    parser.add_argument(
        "--precision",
        action="store",
        default=7,
        type=int,
        help=("Define the number of decimal used in comparing" "float."),
        dest="precision",
    )

    return vars(parser.parse_args())


def _build_test_index(root_path, tests):
    """Look for a valid xml in all the specified directories.

    Check if there are valid xml in all the specified folders. In case append
    the test to a list. If a folder contains more than a single xml, that folder
    will be skipped and a warning is printed. If a folder do not contain any xml
    that folder will be skipped without any message.
    Once all the folder are examinated, a list of all the found tests is printed
    at stdout.

    Print to stdout a list of all the tests found.

    Args:
        root_path: Parent folder of all the tests.
        tests: A sequence of subfolder of root_path to filter the number of
            tests found.

    Returns:
        test_list: A sequence containing all the valid tests found. Each
            element is a tuple of the form (<test_name>, <abs_path_to_xml>)

    """

    root_path = os.path.abspath(root_path)
    tests_folder = [os.path.join(root_path, test) for test in tests]
    test_list = []
    msg = "**Tests that will be executed:**\n"

    for test in tests_folder:
        if not os.path.exists(test) or not os.path.isdir(test):
            raise RuntimeError("Folder %s does not exist!" % test)
        for root, junk, files in os.walk(test):  # pylint: disable=unused-variable
            test_name = os.path.relpath(root, root_path)
            xml_files = [x for x in files if x.endswith(".xml")]
            if len(xml_files) > 1:
                sys.stderr.write(
                    "!W! Skipping test %s: too many xml files!\n" % test_name
                )
                continue
            elif len(xml_files) < 1:
                continue

            if _file_is_test(os.path.join(root, xml_files[0])):
                test_list.append((test_name, os.path.join(root, xml_files[0])))
                msg += " > " + str(os.path.split(root)[1]) + "\n"

    if len(test_list) < 1:
        print("**No test found!**")
    else:
        print(msg)

    return test_list


def _file_is_test(path_to_test):
    """Check if an xml file can be used as regtest.

    To be a valid regtest the xml file must contain the line to specify which
    is the command to be executed as a driver.
    """

    with open(path_to_test) as _file:
        _text = _file.read()
    print(_text[:100])
    return len([x.group(1) for x in REGTEST_STRING_RGX.finditer(_text)]) > 0


class Test(threading.Thread):
    """Contains all the methods used to create, run and compare a test.

    Args:
        index: An integer used to ensure no-overlap between sockets.
        path: The path to the reference xml file.
        name: The path of the test relative to the test_folder.
        root_run: The path where the test will be run.
    """

    def __init__(self, *args, **kwds):
        threading.Thread.__init__(self)
        # multiprocessing.Process.__init__(self)

        self.save_args = kwds

        self.index = kwds["index"]
        self.test_path = kwds["path"]
        self.name = kwds["name"]
        self.run_path = os.path.abspath(kwds["root_run"])

        self.ref_path = os.path.dirname(self.test_path)
        self.filename = os.path.basename(self.test_path)
        self.msg = ""
        self.driver_command = None
        self.input_dir = None
        self.io_dir = None
        self._test_status = "PASSED"
        self.needed_files = None

        self.generate_output = False
        self.compare_output = False
        self.copy_reference = False
        self.die = False

    def copy(self):
        """Useful to create a copy of the actual state of the class.

        Returns:
            A complete copy of the class in his actual state.
        """
        _copy = Test(**self.save_args)
        _copy.ref_path = self.ref_path
        _copy.filename = self.filename
        _copy.msg = self.msg
        _copy.driver_command = self.driver_command
        _copy.input_dir = self.input_dir
        _copy.io_dir = self.io_dir
        _copy.test_status = self.test_status
        _copy.compare_output = self.compare_output
        _copy.generate_output = self.generate_output
        return _copy

    @property
    def test_status(self):
        """Status of the test: always print the most important status!

        An error has always the precedence, then a failure and, only if nothing
        bad happen, print PASSED.
        """
        return self._test_status

    @test_status.setter
    def test_status(self, status):
        status_priority = {"ERROR": 10, "FAILED": 5, "COPIED": 2, "PASSED": 1}

        if status_priority[status] > status_priority[self._test_status]:
            self._test_status = status

    def init_env(self):
        """Prepare the test run.

        Do the following thing:
        * Create the necessary folders.
        * Retrieve the command to run the driver from the xml.
        * Copy all the reference files into the 'input' folder.
        * Copy all the necessary files to the 'io' folder.
        * Change the address of the socket to be unique

        There is some rule on the way the socket address is changed and it is
        important to know this rules when preaparing the test:
        * If unix socket then the <address> tag is changed appending a number
          defined using the self.index.
        * If inet socket then the <port> tag to a number defined using the
          self.index.
        * To ensure that the driver is using the right address the first word
          of the driver command that match the address of the xml file is
          replaced with the new address.
        """

        # Create the necessary folders
        path = self.run_path
        for folder in self.name.split("/"):
            path = os.path.join(path, folder)
            create_dir(path, ignore=True)
        self.run_path = path

        self.io_dir = os.path.join(self.run_path, "io")
        self.input_dir = os.path.join(self.run_path, "input")
        create_dir(self.io_dir, ignore=True)

        # Retrieve the command to run the driver and the needed files
        self.driver_command, self.needed_files = parse_regtest_string(self.test_path)

        # Copy all the reference files to the 'input' folder
        shutil.copytree(self.ref_path, self.input_dir)

        # Copy all only the needed files to the 'io' folder
        for _buffer in self.needed_files:
            try:
                shutil.copy2(os.path.join(self.input_dir, _buffer), self.io_dir)
            except IOError as _err:
                if _err.errno == os.errno.ENOENT:
                    self.test_status = "ERROR"
                    self.msg += "Needed file %s not found!\n" % _buffer

        # Make sure to do not change the reference
        self.test_path = os.path.join(self.io_dir, self.filename)

        # Change the address of the socket to be unique
        xml = etree.parse(self.test_path)
        address_index = self.index
        for ffsocket in xml.findall("./ffsocket"):
            socket_mode = ffsocket.attrib["mode"]
            if socket_mode.lower() == "unix":
                address = ffsocket.find("./address").text.strip()
                new_address = " %s%i " % (address, address_index)
                ffsocket.find("./address").text = new_address
            else:
                address = ffsocket.find("./port").text.strip()
                new_address = "%i" % address_index
                ffsocket.find("./port").text = new_address

            # Change the address used by the driver too!
            for _ii, _buffer in enumerate(self.driver_command):
                # Determine if the address is in a file
                for _word in _buffer.split():
                    _word = os.path.join(self.io_dir, _word)
                    if os.path.exists(_word) and os.path.isfile(_word):
                        inplace_change(_word, address, new_address)

                # Replace only the first occurrence of 'address' in the
                # +driver command!
                self.driver_command[_ii] = _buffer.replace(address, new_address, 1)

            # Since there could be more than a single socket used within a
            # +single simulation, it is important to have a different address for
            # +each socket.
            address_index += 1

        xml.write(self.test_path)

    def run(self):
        if self.generate_output:
            self.init_env()
            # Change to the io directory...
            os.chdir(self.io_dir)
            self._run_ipi()
        elif self.copy_reference:
            # Change to the io directory...
            os.chdir(self.io_dir)
            self.create_reference()
            self.print_report()
        elif self.compare_output:
            # Change to the io directory...
            os.chdir(self.io_dir)
            self._compare()
            self.print_report()

    def _run_ipi(self):
        # Move to the right directory

        driver_out_path = os.path.join(self.io_dir, "driver_output.out")

        timeout_driver = _parser()["maxtime"]
        timeout_ipi = TIMEOUT_IPI
        ipi_command = "i-pi"

        # Run the i-pi code
        ipi_command = shlex.split(ipi_command + " " + self.test_path)
        with open(os.path.join(self.io_dir, "ipi_output.out"), "w") as ipi_out:
            ipi_out.write("*REGTEST* IPI COMMAND: %s\n" % " ".join(ipi_command))
            try:
                ipi_proc = sbps.Popen(
                    ipi_command, bufsize=0, stdout=ipi_out, stderr=sbps.PIPE
                )
            except OSError as _err:
                if _err.errno == os.errno.ENOENT:
                    self.msg += "i-pi command not found!\n"
                    self.test_status = "ERROR"

        # Sleep few seconds waiting for the ipi server start
        time.sleep(IPI_WAITING_TIME)

        # Run the driver code
        driver_prcs = []
        oldpwd = os.getcwd()
        for cmd in self.driver_command:
            instance_folder = os.path.join(self.io_dir, "driver-%i" % len(driver_prcs))
            create_dir(instance_folder, ignore=True)
            for _file in self.needed_files:
                _src = os.path.join(self.io_dir, _file)
                shutil.copy2(_src, instance_folder)

            os.chdir(instance_folder)
            cmd = shlex.split(cmd)
            with open(driver_out_path, "w") as driver_out:
                driver_out.write("*REGTEST* DRIVER COMMAND: %s\n" % " ".join(cmd))
                try:
                    driver_prcs.append(
                        sbps.Popen(
                            cmd,
                            bufsize=0,
                            stdin=None,
                            stdout=driver_out,
                            stderr=sbps.STDOUT,
                        )
                    )
                except OSError as _err:
                    if _err.errno == os.errno.ENOENT:
                        self.msg += "driver command %s not found!\n" % " ".join(cmd)
                        self.test_status = "ERROR"
            os.chdir(oldpwd)

        init_time = time.time()
        finished = 0
        while finished < len(driver_prcs):
            if self.die:
                timeout_driver = -1000
            for prc in driver_prcs:
                if prc.poll() is not None:
                    finished += 1
            time.sleep(0.5)
            time_to_stop = timeout_driver + init_time - time.time()
            if time_to_stop < -0.5:
                for prc in driver_prcs:
                    if prc.poll() is None:
                        prc.terminate()
                    finished += 1
                self.test_status = "ERROR"
                self.msg += "The drivers took too long:\n {:d}s > {:d}s\n".format(
                    int(TIMEOUT_DRIVER - time_to_stop), int(TIMEOUT_DRIVER)
                )

        init_time = time.time()
        while ipi_proc.poll() is None:
            if self.die:
                timeout_ipi = -1000
            time.sleep(0.5)
            time_to_stop = timeout_ipi + init_time - time.time()
            if time_to_stop < -0.5:
                ipi_proc.terminate()
                self.test_status = "ERROR"
                self.msg += "i-PI took too long after the driver finished!\n"
                time.sleep(11)  # i-pi took approximately 10 sec to exit

        stdout, stderr = ipi_proc.communicate()  # pylint: disable=unused-variable
        if len(stderr) != 0:
            self.test_status = "ERROR"
            self.msg += "\ni-PI produced the following error:\n"
            self.msg += stderr
            self.msg += "\n"

        return

    def create_reference(self):
        """Determines and passes the reference files to a reference folder."""
        if self.test_status == "PASSED":
            reference_dir = os.path.join(self.ref_path, "regtest-ref")
            create_dir(reference_dir)
            ltraj, lprop = self.get_filesname(
                self.test_path, reference_dir, self.io_dir
            )

            try:
                for prop in lprop:
                    for sprop in prop:
                        remove_file(sprop["old_filename"])
                        shutil.copy2(sprop["new_filename"], sprop["old_filename"])
                for traj in ltraj:
                    for straj in traj:
                        remove_file(straj["old_filename"])
                        shutil.copy2(straj["new_filename"], straj["old_filename"])
            except IOError as e:
                self.test_status = "ERROR"
                self.msg += "Error while copying the new reference!!\n"
                self.msg += "Unable to copy file. %s" % e
            else:
                self.test_status = "COPIED"
        else:
            self.msg += "Errors occured: using this run as reference " "is not safe!\n"

    def print_report(self):
        """This function prints information about the test on stdout.

        The test can have 3 results:
        - PASSED: The computed values match the reference.
        - FAILED: The computed values DO NOT match the reference.
        - ERROR: The test did not provide any results because some errors
            raised while running i-pi or while running the driver.
        There is also another possible result:
        - COPIED: If the user is creating the reference folder, when
            everything goes fine the reported result is "COPIED".

        When something goes wrong information about what is the problem are
        also printed.
        """

        if self.test_status == "ERROR":
            _format = RESET + "%-30s --> " + RED + "%15s" + INFO + " Info: %s" + RESET
        elif self.test_status == "FAILED":
            _format = (
                RESET + "-%30s --> " + YELLOW + "%15s" + INFO + " Info: %s" + RESET
            )
        elif self.test_status == "PASSED":
            _format = RESET + "-%30s --> " + GREEN + "%15s" + RESET
        elif self.test_status == "COPIED":
            _format = RESET + "-%30s --> " + GREEN + "%15s" + RESET

        if len(self.msg) > 0:
            lines = self.msg.split("\n")
            for _ii, _buffer in enumerate(lines):
                if _ii == 0:
                    continue
                lines[_ii] = " " * 57 + _buffer
            _msg = "\n".join(lines)
            msg = _format % (self.name, self.test_status, _msg)
        else:
            _format = "%30s -->  %15s\n"
            msg = _format % (self.name, self.test_status)
        sys.stdout.write(msg)
        # print msg

    def compare_property_files(self, lprop):
        """Comparing property files.

        The files are loaded with np.loadtxt and all the intestation are
        ignored: only the actual values are compared. Thus, the ipi input must
        define the same columns in the same order.
        """

        for prop in lprop:
            for sprop in prop:
                try:
                    old_content = np.loadtxt(sprop["old_filename"])
                except IOError as _err:
                    if _err.errno == os.errno.ENOENT:
                        self.msg += "File %s not found!\n" % sprop["old_filename"]
                        self.test_status = "ERROR"
                        continue

                try:
                    new_content = np.loadtxt(sprop["new_filename"])
                except IOError as _err:
                    if _err.errno == os.errno.ENOENT:
                        self.msg += "File %s not found!\n" % sprop["old_filename"]
                        self.test_status = "ERROR"
                        continue

                try:
                    npt.assert_array_almost_equal(
                        old_content, new_content, _parser()["precision"]
                    )
                except AssertionError:
                    name = os.path.basename(sprop["old_filename"])
                    self.msg += "Differences in the %s file\n" % name
                    self.test_status = "FAILED"
                    continue

    def compare_trajectory_files(self, ltraj):
        """Function to compare trajectory files.

        The idea is to store all the numbers in the file in a list and then using
        numpy to compare the two lists. The numbers are recognized exploiting the
        float function error when applied on strings.

        The strings are compared directly.
        """

        err = False
        for traj in ltraj:
            for straj in traj:
                new_w_list = []
                old_w_list = []
                name = os.path.basename(straj["old_filename"])
                try:
                    old_content = open(straj["old_filename"])
                except IOError as _err:
                    if _err.errno == os.errno.ENOENT:
                        self.msg += "File %s not found!\n" % straj["old_filename"]
                        self.test_status = "ERROR"
                        continue

                try:
                    new_content = open(straj["new_filename"])
                except IOError as _err:
                    if _err.errno == os.errno.ENOENT:
                        self.msg += "File %s not found!\n" % straj["new_filename"]
                        self.test_status = "ERROR"
                        continue

                line_c = 1
                for old_line, new_line in zip(old_content, new_content):
                    word_c = 1
                    for old_w, new_w in zip(old_line.split(), new_line.split()):
                        try:
                            old_w_list.append(float(old_w))
                            new_w_list.append(float(new_w))
                        except ValueError:
                            try:
                                assert old_w == new_w
                            except AssertionError:
                                self.msg += (
                                    "Differences at line %d word %d of file  %s"
                                    % (line_c, word_c, name)
                                )
                                self.test_status = "FAILED"

                        word_c += 1
                    line_c += 1

                try:
                    npt.assert_array_almost_equal(
                        np.array(new_w_list),
                        np.array(old_w_list),
                        _parser()["precision"],
                    )
                except AssertionError:
                    self.msg += "Differences in the %s file\n" % name
                    self.test_status = "FAILED"
                    continue

        return err

    def _compare(self):
        """This is the function that compares all the ipi output.

        The name of the files to compare come from ltraj and lprop.
        """
        ltraj, lprop = self.get_filesname(
            self.test_path, os.path.join(self.input_dir, "regtest-ref"), self.io_dir
        )

        self.compare_property_files(lprop)

        self.compare_trajectory_files(ltraj)

    def get_filesname(self, xml_path, olddir, newdir):
        """The test results should be analyzed number by numbers.

        The idea is that the testing input should never change, then the files
        should be always the same. It would probably be better, anyway, to use
        the i-pi infrastructure to retrieve the right position of the data. In
        fact, this would work as a further testing.

        Args:
            olddir: The path used for the 'old_filename' in the returned
                dictionary.
            newdir: The path used for the 'new_filename' in the returned
                dictionary.

        Returns:
            lprop
            nprop
        """
        # Avoid to print i-pi output
        devnull = open("/dev/null", "w")
        oldstdout_fno = os.dup(sys.stdout.fileno())
        os.dup2(devnull.fileno(), 1)

        # opens & parses the input file

        # get in the input file location so it can find other input files for initialization
        cwd = os.getcwd()
        iodir = os.path.dirname(os.path.realpath(xml_path))
        os.chdir(iodir)

        # print "READING FILE FROM ", iodir
        # print " WHILE RUNNING IN ", cwd
        # print "I have changed directory to ", os.getcwd()

        ifile = open(xml_path, "r")
        xmlrestart = io_xml.xml_parse_file(ifile)  # Parses the file.
        ifile.close()

        isimul = InputSimulation()
        isimul.parse(xmlrestart.fields[0][1])

        simul = isimul.fetch()
        os.chdir(cwd)

        # reconstructs the list of the property and trajectory files
        lprop = []  # list of property files
        ltraj = []  # list of trajectory files
        for o in simul.outtemplate:
            # properties and trajectories are output per system
            if isinstance(o, CheckpointOutput):
                pass
            elif isinstance(o, PropertyOutput):
                nprop = []
                isys = 0
                for _ss in simul.syslist:  # create multiple copies
                    filename = _ss.prefix + o.filename
                    nprop.append(
                        {
                            "old_filename": os.path.join(olddir, filename),
                            "new_filename": os.path.join(newdir, filename),
                            "stride": o.stride,
                            "properties": o.outlist,
                        }
                    )
                    isys += 1
                lprop.append(nprop)

            # trajectories are more complex, as some have per-bead output
            elif isinstance(o, TrajectoryOutput):
                if getkey(o.what) in [
                    "positions",
                    "velocities",
                    "forces",
                    "forces_sc",
                    "extras",
                ]:  # multiple beads
                    nbeads = simul.syslist[0].beads.nbeads
                    for _bi in range(nbeads):
                        ntraj = []
                        isys = 0
                        # zero-padded bead number
                        padb = (
                            "%0"
                            + str(int(1 + np.floor(np.log(nbeads) / np.log(10))))
                            + "d"
                        ) % (_bi)

                        for _ss in simul.syslist:
                            if o.ibead < 0 or o.ibead == _bi:
                                if getkey(o.what) == "extras":
                                    filename = _ss.prefix + o.filename + "_" + padb
                                else:
                                    filename = (
                                        _ss.prefix
                                        + o.filename
                                        + "_"
                                        + padb
                                        + "."
                                        + o.format
                                    )
                                ntraj.append(
                                    {
                                        "old_filename": os.path.join(olddir, filename),
                                        "format": o.format,
                                        "new_filename": os.path.join(newdir, filename),
                                        "stride": o.stride,
                                        "what": o.what,
                                    }
                                )
                            isys += 1
                        if ntraj != []:
                            ltraj.append(ntraj)

                else:
                    ntraj = []
                    isys = 0
                    for _ss in simul.syslist:  # create multiple copies
                        filename = _ss.prefix + o.filename
                        filename = filename + "." + o.format
                        ntraj.append(
                            {
                                "old_filename": os.path.join(olddir, filename),
                                "new_filename": os.path.join(newdir, filename),
                                "format": o.format,
                                "stride": o.stride,
                            }
                        )

                        isys += 1
                    ltraj.append(ntraj)

        os.dup2(oldstdout_fno, 1)
        return ltraj, lprop


#
# Tools Functions ###
#


def remove_file(path):
    """Remove a path only if it exists and is a file!"""
    if os.path.exists(path) and os.path.isfile(path):
        os.remove(path)


def create_dir(folder_path, ignore=False):
    """Create a folder after user consense.

    If the path asked by the user is not existing, just create that folder
    otherwise asks the user his opinion on deleting the existing folder/file.

    Args:
        folder_path: The path of the folder that should be created.
        ignore: If True and the path already exists, go on...

    """
    if os.path.exists(folder_path):
        if ignore:
            return True
        if answer_is_y(
            "!W! %s already exists!\n    Do you want to delete it and "
            "proceed?[y/n]\n" % folder_path
        ):
            try:
                if os.path.isdir(folder_path):
                    shutil.rmtree(folder_path)
                elif os.path.isfile(folder_path):
                    os.remove(folder_path)
                else:
                    raise RuntimeError
            except OSError:
                raise RuntimeError(
                    "I cannot remove the file. " "Try manually and restart this script!"
                )
        else:
            raise SystemExit("User rules!")

    os.mkdir(folder_path)

    return True


def answer_is_y(msg):
    """A simple function to interrogate the user on y/n questions.

    Only 'yes' and 'y' are counted as yes answer and only 'n' and 'no' are
    valied negative answer. All the answer are case insensitive. The function
    will ask for an answer until the user do not reply with a valid character.

    Args:
        msg: The question...

    Return:
        True if the user answer yes or False if the user answer no.
    """

    _yes = ["yes", "y"]
    _no = ["no", "n"]
    if not msg.endswith("\n"):
        msg += "\n"
    if not msg.startswith("\n"):
        msg = "\n" + msg

    answer = ""
    while answer.lower() not in _yes and answer not in _no:
        sys.stderr.write(msg)

        answer = input()

        if answer.lower() in _yes:
            return True
        elif answer.lower() in _no:
            return False


def parse_regtest_string(test_path):
    """Retrieve the commands and the dependecies from the xml input.

    To run the test, the command of the driver and the necessary files must be
    known. These information are enclosed in a comment into the xml file. This
    function parses this comment and return the command to be use as driver and
    a list of files needed.

    Args:
        test_path: Path to the xml file.

    Returns:
        commands: list of command to be used to run the driver.
        dependencies: list of file needed to run the test (the xml file is
            already obviuous).
    """
    command_string_rgx = re.compile(
        r"^COMMAND\(?(\d*)\)?\s*([ \w.\+\-\(\)\<\>\:]*)$", re.MULTILINE
    )
    dependency_string_rgx = re.compile(
        r"^DEPENDENCIES\s*([ \w.\+\-\(\)]*)$", re.MULTILINE
    )

    with open(test_path) as _buffer:
        _text = _buffer.read()

    regtest_string = []
    dependencies = [test_path]
    commands = []

    regtest_string = REGTEST_STRING_RGX.findall(_text)

    for _xx in dependency_string_rgx.finditer("\n".join(regtest_string)):
        dependencies += _xx.group(1).split()
    for _xx in command_string_rgx.finditer("\n".join(regtest_string)):
        try:
            commands += [_xx.group(2)] * int(_xx.group(1))
        except ValueError:
            commands += [_xx.group(2)]

    return commands, dependencies


def inplace_change(filename, old_string, new_string):
    """Replace a string in a file.

    Replace 'old_string' with 'new_string' into 'filename'.

    Args:
        filename: The filename where the string must be replaced.
        old_string: The string to replace.
        new_string: The string that will be replaced.

    Returns:
        Will be True if something has been replaced, False otherwise.
    """
    # Safely read the input filename using 'with'
    with open(filename) as _file:
        _content = _file.read()
        if old_string not in _content:
            return False
        if new_string in _content:
            return True

    # Safely write the changed content, if found in the file
    with open(filename, "w") as _file:
        _content = _content.replace(old_string, new_string)
        _file.write(_content)
        return True


if __name__ == "__main__":
    main()
