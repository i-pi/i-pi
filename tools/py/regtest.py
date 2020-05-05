#!/usr/bin/env python2

"""
This script runs i-PI regression tests.

regtest by default searches recursively for i-PI regression tests in the
current directory and executes them in ./regtest-run

With option --create-reference it does not run tests but prepares
reference output files for regression testing in the regrest-ref
directory for each test case.

There are also following options are available:
  --test-cases-directory TEST_CASES_DIRECTORY
    regtest will recursively search for regression tests in the given
    directory instead of the default current directory
  --run-directory RUN_DIRECTORY
    regtest will run test instance in the directory RUN_DIRECTORY/regtest-run

Creating a regression test case
===============================
In order to create a regression test case, modify the i-PI input: add a regtest header. Also, do not forget to change the number of steps in the i-PI input.

Input must contain a header that specifies which
files are needed to run the test case (dependencies) and what command should be used to run the
driver. You must follow the exact syntax (only whitespace does not matter).

Here is the example of the header:

<!--REGTEST
DEPENDENCIES  h5o2.dms4B.coeff.com.dat h5o2.pes4B.coeff.dat h5o2+.xyz
COMMAND(10)    i-pi-driver -u -h REGTEST_SOCKET -m zundel
COMMAND       i-pi-driver -u -h REGTEST_SOCKET -m zundel
ENDREGTEST-->

`DEPENDECIES`: files needed to run the test (except the i-PI input)
`COMMAND(n)`: command used to run the driver. The `n` is the number of
    instances of the COMMAND to be created. In the example there will be a
    total of 11 instances of the driver running at once.

This script will replace socket adresses. In case the name of the socket
matches more than once in the COMMAND instruction above, only the first one will be
changed. For example: using the command "i-pi-driver -u -m zundel -h zundel"
would only replace the first string "zundel" and cause an error.
"""

import sys
import re
import time
import os
from copy import deepcopy
import shutil
from collections import deque
import xml.etree.ElementTree as etree
import shlex
import subprocess
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml
from ipi.engine.outputs import CheckpointOutput, PropertyOutput, TrajectoryOutput
from ipi.engine.properties import getkey
import numpy as np
import numpy.testing as npt
import argparse
import traceback


def main():
    arguments = command_parser()
    run_dir = os.path.join(
        os.path.abspath(arguments["run_directory"]), Parameters.run_directory
    )
    test_dir = os.path.abspath(arguments["test_cases_directory"])
    is_in_reference_mode = arguments["is_in_reference_mode"]
    try:
        os.makedirs(run_dir)
    except OSError:
        print(
            "The directory %s exists. Do you want to delete its contents and continue? (y/n)"
            % run_dir
        )
        if answer_is_y():
            shutil.rmtree(run_dir)
            os.makedirs(run_dir)
        else:
            quit("Terminated")
    if is_in_reference_mode:
        print("Do you agree to replace references if they exist? (y/n)")
        if answer_is_y():
            replace_references = True
        else:
            replace_references = False

    list_of_test_candidates = find_test_candidates(test_dir)
    list_of_test_cases = []
    for test_candidate in list_of_test_candidates:
        try:
            list_of_test_cases.append(TestCase(test_candidate))
        except TypeError as e:
            print("--> Could not create test case:", test_candidate.name)
            print(str(e))
            continue
    queue_of_test_instances = deque([])
    counter = Counter()
    print()
    print("List of test cases to be executed:")
    for test_case in list_of_test_cases:
        test_instance_dir = os.path.join(run_dir, test_case.name)
        queue_of_test_instances.append(
            TestInstance(test_case, test_instance_dir, counter)
        )
        print("-->", test_case.name)
    print()
    print("Executing following test instances:")
    while queue_of_test_instances:
        test_instance = queue_of_test_instances.popleft()
        try:
            sys.stdout.write("--> " + test_instance.name)
            sys.stdout.flush()
            test_instance.run()
            if not is_in_reference_mode:
                try:
                    test_passed = True
                    differences = test_instance.compare_output_with_reference()
                    for results in differences:
                        if not results.files_are_equal():
                            test_passed = False
                    if test_passed:
                        print("    PASSED")
                    else:
                        print("    FAILED")
                        for results in differences:
                            results.print_differences()
                except ValueError as e:
                    print(str(e), file=sys.stderr)
            else:
                try:
                    test_instance.put_output_as_reference()
                except ValueError:
                    if replace_references:
                        test_instance.test_case.remove_reference()
                        test_instance.put_output_as_reference()
                        print("    SUCCESS: References replaced")
                    else:
                        print("    SKIPPING: References not copied")

        except (RegtestTimeout, WrongDriverCommand, IPIError, OutputCorrupted) as e:
            print("    ERROR")
            print(
                str(e),
                "Test instance run directory:",
                test_instance.run_dir,
                file=sys.stderr,
            )

            continue
        except WrongIPICommand as e:
            sys.exit(str(e))


class TestCase:

    """
    Stores test case, paths to original files which define a regtest.

    Attributes:
        input_path: path to test input file
        name: test case name based on its relative directory
        root: directory where test input is stored
        dependencies_list: list of paths to all dependent files
        dependencies_dir: direcotry where dependent files are stored
        reference_dir: directory where reference output is stored
    """

    def __init__(self, test_candidate):
        """
        Before initializing TestCase based on TestCandidate,
        checks if dependencies exist.
        Arguments:
            test_candidate: instance of TestCandidate
        """
        dependencies_dir = os.path.dirname(test_candidate.input_path)
        commands, dependencies = parse_regtest_string(test_candidate.input_path)
        dependencies_list = [os.path.join(dependencies_dir, x) for x in dependencies]
        if not all([os.path.exists(x) for x in dependencies_list]):
            raise TypeError(
                "Dependencies for file: "
                + str(test_candidate.input_path)
                + " absent in: "
                + str(dependencies_dir)
            )
        self.input_path = test_candidate.input_path
        self.name = test_candidate.name
        self.root = os.path.dirname(test_candidate.input_path)
        self.dependencies_list = dependencies_list
        self.dependencies_dir = dependencies_dir
        self.reference_dir = os.path.join(self.root, Parameters.reference_directory)

    def how_many_sockets(self):
        """
        Checks how many sockets the test case needs.
        Returns:
            sockets: number of sockets needed
        """
        xml = etree.parse(self.input_path)
        sockets = 0
        for ffsocket in xml.findall("./ffsocket"):
            sockets += 1
        return sockets

    def get_reference_output(self):
        """
        Returns:
            TestOutput object initialized from test case's reference.
        """
        return TestOutput(self.input_path, self.reference_dir)

    def remove_reference(self):
        """
        Removes reference belonging to the test case.
        """
        shutil.rmtree(self.reference_dir)


class Counter:

    """
    Counter class for attributing socket numbers.
    """

    socket = 10001

    def attribute_socket(self):
        socket = Counter.socket
        Counter.socket += 1
        return socket


class TestInstance:

    """
    Stores test instance which is build from a test case.

    Attributes:
        test_case: test case on which test instance is based
        path_ipi_input: i-pi input of the test instance
        name: name of the test instance (the same as the test case)
        driver_signatures: tuples of driver commands and their directories
        run_dir: root directory of the test instance
        dependencies: absolute paths of dependencies
        ipi_command: command to run the test instance
    """

    def __init__(self, test_case, dir, counter):
        """
        Before initializing instance from TestCase, checks if the path
        to directory does not exist. Also creates modified copies of
        dependencies with correct socket names. Updaets driver commands.

        Arguments:
            test_case: TestCase object
            dir: root directory of the test instance
            counter: Counter object
        """
        run_dir = os.path.abspath(dir)
        try:
            os.makedirs(run_dir)
        except OSError:
            if os.path.exists(run_dir):
                raise ValueError(
                    "Directory %s exists. The dir parameter must be a path to nonexistent directory"
                    % str(run_dir)
                )
            else:
                raise
        test_dependencies = []
        for files in test_case.dependencies_list:
            shutil.copy(files, run_dir)
            test_dependencies.append(os.path.join(run_dir, os.path.basename(files)))

        xml = etree.parse(test_case.input_path)
        commands, dependencies = parse_regtest_string(test_case.input_path)
        for ffsocket in xml.findall("./ffsocket"):
            address_index = counter.attribute_socket()
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
            for _ii, _buffer in enumerate(commands):
                # Determine if the address is in a file
                for _word in _buffer.split():
                    _word = os.path.join(run_dir, _word)
                    if os.path.exists(_word) and os.path.isfile(_word):
                        inplace_change(_word, address, new_address)

                # Replace only the first occurrence of 'address' in the
                # driver command!
                commands[_ii] = _buffer.replace(address, new_address, 1)

        xml.write(os.path.join(run_dir, "input.xml"))
        driver_signatures = []  # tuple: command + dir
        for cmd in commands:
            instance_folder = os.path.join(
                run_dir, "driver-%i" % len(driver_signatures)
            )
            os.makedirs(instance_folder)
            for _file in test_dependencies:
                shutil.copy(_file, instance_folder)

            driver_command = shlex.split(cmd)
            driver_signatures.append((driver_command, instance_folder))
        self.test_case = test_case
        self.path_ipi_input = os.path.join(run_dir, "input.xml")
        self.name = test_case.name
        self.driver_signatures = driver_signatures
        self.run_dir = run_dir
        self.dependencies = test_dependencies
        self.ipi_command = shlex.split(
            "i-pi" + " " + os.path.join(run_dir, "input.xml")
        )

    def run(self):
        """
        Runs the test instance.
        """
        ipi_output_path = os.path.join(self.run_dir, Parameters.ipi_output_file)
        # Run the i-pi code
        os.chdir(self.run_dir)
        driver_prcs = deque([])
        try:
            with open(ipi_output_path, "w") as ipi_out:
                ipi_out.write(
                    "*REGTEST* IPI COMMAND: %s\n" % " ".join(self.ipi_command)
                )
                try:
                    ipi_proc = subprocess.Popen(
                        self.ipi_command,
                        bufsize=0,
                        stdout=ipi_out,
                        stderr=subprocess.PIPE,
                    )
                except OSError as _err:
                    if _err.errno == os.errno.ENOENT:
                        raise WrongIPICommand(
                            "i-pi command not found!", self.ipi_command
                        )
                    else:
                        raise
            time.sleep(Parameters.ipi_starting_time)

            if ipi_proc.poll() is not None:
                stdout, stderr = ipi_proc.communicate()
                raise IPIError("I-PI error on start\n" + stderr)

            # Run the driver code
            for cmd, instance_folder in self.driver_signatures:
                os.chdir(instance_folder)
                driver_out_path = os.path.join(instance_folder, "driver.out")
                with open(driver_out_path, "w") as driver_out:
                    driver_out.write("*REGTEST* DRIVER COMMAND: %s\n" % " ".join(cmd))
                    try:
                        driver_prcs.append(
                            subprocess.Popen(
                                cmd,
                                bufsize=0,
                                stdin=None,
                                stdout=driver_out,
                                stderr=subprocess.PIPE,
                            )
                        )
                    except OSError as _err:
                        if _err.errno == os.errno.ENOENT:
                            raise WrongDriverCommand(
                                "driver command %s not found!\n" % " ".join(cmd)
                            )
                        else:
                            raise

            driver_init_time = time.time()
            driver_errors = []
            driver_return_codes = []
            # while queue is not empty
            while driver_prcs:
                # take first process in the queue
                prc = driver_prcs.popleft()
                # if driver process is running
                if prc.poll() is None:
                    # enqueue it and wait a bit before taking next one
                    driver_prcs.append(prc)
                    time.sleep(Parameters.sleep_time)
                else:
                    driver_return_codes.append(prc.returncode)
                    stdout, stderr = prc.communicate()
                    if stderr:
                        driver_errors.append(stderr)
                time_elapsed = time.time() - driver_init_time
                if time_elapsed > Parameters.driver_timeout:
                    for prc in driver_prcs:
                        if prc.poll() is None:
                            prc.terminate()
                    ipi_proc.terminate()
                    stdout, stderr = ipi_proc.communicate()
                    if stderr:
                        with open(ipi_output_path, "a") as ipi_out:
                            ipi_out.write(stderr)
                        raise IPIError("I-PI Error\n")
                    raise RegtestTimeout("Driver timeout")

            drivers_terminated_time = time.time()
            while ipi_proc.poll() is None:
                time.sleep(Parameters.sleep_time)
                time_elapsed_after_drivers = time.time() - drivers_terminated_time
                if time_elapsed_after_drivers > Parameters.ipi_shutdown_time:
                    ipi_proc.terminate()
                    stdout, stderr = ipi_proc.communicate()
                    if stderr:
                        with open(ipi_output_path, "a") as ipi_out:
                            ipi_out.write(stderr)
                        raise IPIError("i-PI Error\n" + stderr)
                    raise RegtestTimeout(
                        "i-PI was terminated by regtest after %f s after drivers returned\n"
                        % Parameters.ipi_shutdown_time
                        + "This might be caused by error in driver(s).\n"
                        + "Driver return codes: "
                        + ", ".join(str(x) for x in driver_return_codes)
                        + "\n"
                        "Driver errors:" + "\n".join(driver_errors) + "\n"
                    )

            stdout, stderr = ipi_proc.communicate()
            if stderr:
                with open(ipi_output_path, "a") as ipi_out:
                    ipi_out.write(stderr)
                raise IPIError("I-PI Error\n" + stderr)
        except KeyboardInterrupt as e:
            print("\nKeyboard interrupt!")
            traceback.print_exc(file=sys.stdout)
            if ipi_proc.poll() is None:
                print("Trying to shut down i-PI cleanly...")
                # Let I-PI start before we terminate it, otherwise it can hang up
                time.sleep(Parameters.ipi_starting_time)
                ipi_proc.terminate()
                stdout, stderr = ipi_proc.communicate()
                if stderr:
                    with open(ipi_output_path, "a") as ipi_out:
                        ipi_out.write(stderr)
                print("i-PI terminated")
            for prc in driver_prcs:
                if prc.poll() is None:
                    prc.terminate()
            sys.exit("Regtest termination due to user interruption.")

    def get_output(self):
        """
        Create TestOutput object based on test instance output.
        """
        return TestOutput(self.path_ipi_input, self.run_dir)

    def put_output_as_reference(self):
        """
        Uses output of test instance run and moves the files to
        the reference of the test case.
        """
        output = self.get_output()
        output.put(self.test_case.reference_dir)
        return

    def compare_output_with_reference(self):
        """
        Compares test instance output with test case reference.
        """
        output = self.get_output()
        reference = self.test_case.get_reference_output()
        return output.compare(reference)


class TestOutput:

    """
    Stores paths to test outputs.
    Attributes:
        abs_files: list of absolute paths to test output files
        fileset: set of basenames of test output files
        filetuples: set of tuples (basename, absname) of output files
    """

    def __init__(self, xml, path):
        """
        Args:
            xml: path to test input file
            path: path to test outputs
        """
        abs_xml = os.path.abspath(xml)
        if not os.path.isfile(abs_xml):
            raise ValueError("Input file does not exist %s" % abs_xml)
        filelist = get_output_filenames(abs_xml)
        abs_files = [os.path.join(path, x) for x in filelist]
        filetuples = [(os.path.join(path, x), x) for x in filelist]
        for files in abs_files:
            if not (os.path.isfile(files)):
                raise OutputCorrupted("Expected output file %s not exist" % files)
        self.files = abs_files
        self.fileset = set(filelist)
        self.filetuples = set(filetuples)

    def compare(self, test_output):
        """
        Compares output of test instance with given test_output.
        Arguments:
            test_output: instance of TestOutput
        Returns:
            report: instance of ComparisonResult
        """
        # Check if file lists are the same
        assert self.fileset == test_output.fileset
        # Check if they all exist
        for files in self.files:
            assert os.path.isfile(files)
        for files in test_output.files:
            assert os.path.isfile(files)
        report = []
        for file1 in self.filetuples:
            file2 = [x for x in test_output.filetuples if x[1] == file1[1]]
            differences = compare_files(file1[0], file2[0][0])
            result = ComparisonResult(file1[0], file2[0][0], differences)
            report.append(result)
        return report

    def put(self, path):
        """
        Copies output of test instance to the given path.
        Arguments:
            path: directory to which test outputs will be copied
        """
        try:
            os.makedirs(path)
        except OSError:
            if os.path.exists(path):
                raise ValueError(
                    "Directory %s exists. The dir parameter must be a path to nonexistent directory"
                    % str(path)
                )
            else:
                raise
        for file in self.files:
            shutil.copy(file, path)


class ComparisonResult:

    """
    Stores result of comparing two output files.
    Attributes:
        file1: absolute path of 1st file compared
        file2: absolute path of 2nd file compared
        differences: list of tuples (line, word) specifying differences between file1 and file2
    """

    def __init__(self, file1, file2, differences):
        """
        Initializes a ComparisonResult using arguments
            Arguments:
        file1: absolute path of 1st file compared
        file2: absolute path of 2nd file compared
        differences: list of tuples (line, word) specifying differences between file1 and file2
        """
        self.file1 = file1
        self.file2 = file2
        self.list_of_differences = differences

    def append(self, place):
        """
        Appends a tuple to list of differences
        Arguments:
            place: tuple (line, word)
        """
        self.list_of_differences.append(place)

    def print_differences(self):
        """
        Prints differences between file1 and file2
        """
        print("Files:", self.file1, self.file2)
        for element in self.list_of_differences:
            print("Difference in line", element[0], "in word", element[1])

    def files_are_equal(self):
        """
        Returns True if files are equal.
        """
        return len(self.list_of_differences) == 0


class TestCandidate:

    """
    Contains valid candidates for TestCase.
    Arguments:
        input_path: absolute path to test input
        name: name of the test candidate
    """

    def __init__(self, name, input_path):
        """
        Args:
            name: name of the test candidate
            input_path: absolute path to test input
        """
        input_abspath = os.path.abspath(input_path)
        if file_is_test(os.path.join(input_abspath)):
            self.input_path = input_abspath
            self.name = name
        else:
            raise ValueError("Not a valid xml test")


class Parameters:

    """
    Set of global parameters
    """

    regtest_string = r"<!--\s*REGTEST\s+([\s\w\W]*)\s+ENDREGTEST\s*-->"
    command_string = r"^COMMAND\(?(\d*)\)?\s*([ \w.\+\-\(\)\<\>\:]*)$"
    dependencies_string = r"^DEPENDENCIES\s*([ \w.\+\-\(\)]*)$"
    default_root_run_directory = "."
    run_directory = "regtest-run"
    test_cases_directory = "."
    reference_directory = "regtest-ref"
    precision = 4
    driver_timeout = 600
    ipi_output_file = "ipi_output.out"
    ipi_shutdown_time = 11
    ipi_starting_time = 5
    sleep_time = 0.1

    @staticmethod
    def compile_regtest_string():
        return re.compile(Parameters.regtest_string)

    @staticmethod
    def compile_command_string():
        return re.compile(Parameters.command_string, re.MULTILINE)

    @staticmethod
    def compile_dependencies_string():
        return re.compile(Parameters.dependencies_string, re.MULTILINE)


class RegtestTimeout(Exception):

    """Raise when subprocess of regttest timeout"""


class WrongIPICommand(Exception):

    """Raise when i-pi command is not found"""


class WrongDriverCommand(Exception):

    """Raise when driver command is not found"""


class OutputCorrupted(Exception):

    """Raise when test output (either reference or test) has missing files"""


class IPIError(Exception):

    """Raise when I-PI terminates with error"""


def command_parser():
    """
    Parse the argument lists given as input.

    Returns:
        A dictionary containing all the input options and argument.
    """
    parser = argparse.ArgumentParser(
        description="Recursively looks for "
        "i-pi regression tests in the current directory, "
        "creates './regtest-run' directory and runs tests there. "
        "With the --create-reference flag creates the reference outputs."
    )
    parser.add_argument(
        "--test-cases-directory",
        action="store",
        type=str,
        default=Parameters.test_cases_directory,
        help="The directory which is recursively searched for regression tests",
        dest="test_cases_directory",
    )
    parser.add_argument(
        "--run-directory",
        action="store",
        type=str,
        default=Parameters.default_root_run_directory,
        help=("Directory where the 'regtest-run' " "directory will be created."),
        dest="run_directory",
    )
    parser.add_argument(
        "--create-reference",
        action="store_true",
        default=False,
        help=(
            "Do not compare output and create the relevant "
            "reference outputs for test cases."
        ),
        dest="is_in_reference_mode",
    )

    return vars(parser.parse_args())


def file_is_test(path_of_xml_file):
    """
    Check if an xml file can be used as regtest.
    Args:
        path_of_xml_file: absolute path to regtest file
    Returns:
        True if file is a valid regtest input.
    """

    with open(path_of_xml_file) as _file:
        _text = _file.read()
    regtest_string_rgx = Parameters.compile_regtest_string()
    return len([x.group(1) for x in regtest_string_rgx.finditer(_text)]) > 0


def parse_regtest_string(xml_input_path):
    """
    Retrieve the commands and the dependecies from the xml input.

    This function parses the comment at the beginning of xml file
    and returns the commands to be used as driver and
    a list of dependencies (files required).

    Args:
        xml_input_path: Path of the xml file.

    Returns:
        commands: list of command to be used to run the driver.
        dependencies: list of files needed to run the test
    """

    regtest_string_rgx = Parameters.compile_regtest_string()
    command_string_rgx = Parameters.compile_command_string()
    dependency_string_rgx = Parameters.compile_dependencies_string()

    with open(xml_input_path) as _buffer:
        _text = _buffer.read()

    dependencies = []
    commands = []

    regtest_string = regtest_string_rgx.findall(_text)

    for _xx in dependency_string_rgx.finditer("\n".join(regtest_string)):
        dependencies += _xx.group(1).split()
    for _xx in command_string_rgx.finditer("\n".join(regtest_string)):
        try:
            commands += [_xx.group(2)] * int(_xx.group(1))
        except ValueError:
            commands += [_xx.group(2)]

    return commands, dependencies


def find_test_candidates(root_path):
    """
    Look for a valid xml regtests recursively in the root_path.

    If a directory contains more than a single xml,
    it is skipped and a warning is printed.
    If a directory does not contain any xml, it is skipped without message.

    Args:
        root_path: Parent folder of all the tests.

    Returns:
        list of TestCandidates
    """

    abs_root_path = os.path.abspath(root_path)
    test_list = []

    if not os.path.exists(abs_root_path) or not os.path.isdir(abs_root_path):
        raise RuntimeError("Folder %s does not exist!" % abs_root_path)

    for root, dirs, files in os.walk(abs_root_path):
        test_name = os.path.relpath(root, root_path)
        if test_name == ".":
            test_name = "_root_test_case_"
        xml_files = [x for x in files if x.endswith(".xml")]

        if len(xml_files) > 1:
            # if there are more xml files, it is not clear which one to use
            continue
        elif len(xml_files) == 1:
            try:
                test_list.append(
                    TestCandidate(test_name, os.path.join(root, xml_files[0]))
                )
            except ValueError:
                continue
    return test_list


def check_presence_of_dependencies(xml_file, path):
    """
    Based on regtest input, checks if dependencies are present in path.

    Args:
        xml_file: absolute path to regtest input
        path: path were dependencies should be present

    Returns:
        True if the dependencies are present, False otherwise.
    """
    commands, dependencies = parse_regtest_string(xml_file)
    are_present = True
    abs_dependencies = [os.path.join(os.path.abspath(path), x) for x in dependencies]

    for files in dependencies:
        file_path = os.path.join(path, files)
        if not os.path.isfile(file_path):
            are_present = False
    return are_present


def inplace_change(filename, old_string, new_string):
    """
    Replace 'old_string' with 'new_string' in 'filename'.

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


def get_output_filenames(xml_path):
    """
    Predicts ipi output filenames based on I-PI input file.
    Args:
        xml_path: absolute path to I-PI input file
    Returns:
        list of I-PI output files that will be generated
        running I-PI using this input
    TODO A more elegant solution. This function launches I-PI
    in the directory of the input file and dumps output to dev/null,
    which would make potential bugs hard to detect.
    """
    # Avoid to print i-pi output
    devnull = open("/dev/null", "w")
    oldstdout_fno = os.dup(sys.stdout.fileno())
    os.dup2(devnull.fileno(), 1)

    xml_path = os.path.abspath(xml_path)
    os.chdir(os.path.dirname(xml_path))
    # i-pi xml file parser
    ifile = open(xml_path, "r")
    xmlrestart = io_xml.xml_parse_file(ifile)
    ifile.close()

    isimul = InputSimulation()
    isimul.parse(xmlrestart.fields[0][1])
    simul = isimul.fetch()

    # reconstructs the list of the property and trajectory files
    lprop = []  # list of property files
    ltraj = []  # list of trajectory files
    for o in simul.outtemplate:
        o = deepcopy(o)  # avoids overwriting the actual filename
        if simul.outtemplate.prefix != "":
            o.filename = simul.outtemplate.prefix + "." + o.filename
        # properties and trajectories are output per system
        if isinstance(o, CheckpointOutput):
            pass
        elif isinstance(o, PropertyOutput):
            for _ss in simul.syslist:  # create multiple copies
                filename = o.filename
                if _ss.prefix != "":
                    filename = _ss.prefix + "_" + filename
                lprop.append(filename)
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
                    # zero-padded bead number
                    padb = (
                        "%0" + str(int(1 + np.floor(np.log(nbeads) / np.log(10)))) + "d"
                    ) % (_bi)

                    for _ss in simul.syslist:
                        if (
                            o.ibead < 0 and ((_bi % (-o.ibead) == 0))
                        ) or o.ibead == _bi:
                            filename = o.filename
                            if _ss.prefix != "":
                                filename = _ss.prefix + "_" + filename
                            if getkey(o.what) == "extras":
                                filename += "_" + padb
                            else:
                                filename += "_" + padb + "." + o.format
                            ltraj.append(filename)
            else:
                for _ss in simul.syslist:  # create multiple copies
                    filename = o.filename
                    if _ss.prefix != "":
                        filename = _ss.prefix + "_" + filename
                    filename += "." + o.format
                    ltraj.append(filename)
    os.dup2(oldstdout_fno, 1)
    return ltraj + lprop


def compare_files(file1, file2):
    """
    Compare file1 with file2 and report differences.

    The strings are compared for equality and floats are compared
    to some precision using numpy.

    Args:
        file1, file2: absolute paths to text files

    Return:
        List of tuples of differences (line, word)
    """

    differences = []
    with open(file1, "r") as content_file1:
        content1 = content_file1.readlines()
    with open(file2, "r") as content_file2:
        content2 = content_file2.readlines()

    line_count = 1
    for line_in_file1, line_in_file2 in zip(content1, content2):
        word_count = 1
        for word_in_file1, word_in_file2 in zip(
            line_in_file1.split(), line_in_file2.split()
        ):
            try:
                float_in_file1 = float(word_in_file1)
                float_in_file2 = float(word_in_file2)
                if not np.isclose(float_in_file1, float_in_file2, rtol=1e-4):
                    differences.append((line_count, word_count))
            except ValueError:
                if not word_in_file1 == word_in_file2:
                    differences.append((line_count, word_count))
            word_count += 1
        line_count += 1
    return differences


def answer_is_y():
    """
    Get y/n answer from standard input.

    Only 'yes' and 'y' are counted as yes answer and only 'n' and 'no' are
    valid negative answer. All the answer are case insensitive. The function
    will ask for an answer until the user do not reply with a valid character.

    Return:
        True if the user answer yes or False if the user answer no.
    """

    _yes = ["yes", "y"]
    _no = ["no", "n"]
    answer = input()
    if answer.lower() in _yes:
        return True
    elif answer.lower() in _no:
        return False


if __name__ == "__main__":
    main()
