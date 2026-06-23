import subprocess as sp
import contextlib
import io
import os
from pathlib import Path
import shutil
import signal
import traceback
import xml.etree.ElementTree as ET
import tempfile
import time
import glob

from ipi.scripting import InteractiveSimulation
from ipi.utils.softexit import softexit


def copy_tree(src, dst):  # emulates distutils copy_tree
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


clean_all = False
debug = False
TIMEOUT = 20

fortran_driver_models = [
    "dummy",
    "lj",
    "sg",
    "harm",
    "harm3d",
    "morse",
    "zundel",
    "qtip4pf",
    "qtip4pf-sr",
    "qtip4pf-c-json",
    "qtip4pf-c-1",
    "qtip4pf-c-2",
    "pswater",
    "eckart",
    "ch4hcbe",
    "MB",
    "ljpolymer",
    "doublewell_1D",
    "doublewell",
    "gas",
    "noo3-h2o",
    "water_dip_pol",
]

# We should do this automatically but for now we do it explicitly here
python_driver_models = [
    "dummy",
    "harmonic",
    "DW_friction",
    "DW_explicit",
    "rascal",
    "DoubleWell",
]


def get_test_list(parent, skip_noauto=True):
    """This function recursively searches for test"""
    folders = [x[0] for x in os.walk(parent)]
    reg_tests = list()

    for ff in folders:
        if skip_noauto and ff[-7:].lower() == ".noauto":
            continue
        if os.path.isfile(Path(ff) / "input.xml"):
            reg_tests.append(ff)

    return reg_tests


def clean_tmp_dir():
    if clean_all:
        try:
            for filename in glob.glob("/tmp/ipi_*"):
                os.remove(filename)
        except:
            pass


class InProcessIPI:
    """Popen-like wrapper that runs i-PI through InteractiveSimulation.

    The regression-test runner uses a subprocess-like interface to start i-PI,
    wait for driver sockets and collect output. This wrapper keeps that
    interface while running i-PI in the pytest process, so coverage can observe
    the simulation code.
    """

    def __init__(self, command, cwd):
        """Initializes an interactive simulation from a command.

        Args:
            command: i-PI command line. The last token is interpreted as the
                input file name, matching the current regression-test calls.
            cwd: Directory in which the simulation should be initialized.
        """

        self.command = command
        self.cwd = cwd
        self.stdout = io.StringIO()
        self.stderr = io.StringIO()
        self.returncode = None
        self.simulation = None
        self._cleanup_done = False
        self._ran = False
        old_cwd = os.getcwd()
        try:
            os.chdir(self.cwd)
            with contextlib.redirect_stdout(self.stdout), contextlib.redirect_stderr(
                self.stderr
            ):
                with open(self.command.split()[-1]) as input_file:
                    self.simulation = InteractiveSimulation(input_file)
        except BaseException:
            self.returncode = 1
            self.stderr.write(traceback.format_exc())
            self.kill()
        finally:
            os.chdir(old_cwd)

    def poll(self):
        """Returns the current process-style return code."""

        return self.returncode

    def communicate(self, timeout=None):
        """Runs the simulation and returns captured output.

        Args:
            timeout: Maximum number of seconds to wait before raising
                :class:`subprocess.TimeoutExpired`.

        Returns:
            A tuple ``(stdout, stderr)`` encoded as bytes, matching
            :meth:`subprocess.Popen.communicate`.
        """

        if self._ran or self.returncode is not None:
            self.kill()
            return (
                self.stdout.getvalue().encode("ascii", errors="replace"),
                self.stderr.getvalue().encode("ascii", errors="replace"),
            )

        old_cwd = os.getcwd()
        old_alarm = signal.getsignal(signal.SIGALRM)

        def timeout_handler(signum, frame):
            """Raises the subprocess timeout exception when SIGALRM fires."""

            raise sp.TimeoutExpired(self.command, timeout)

        try:
            os.chdir(self.cwd)
            if timeout is not None:
                signal.signal(signal.SIGALRM, timeout_handler)
                signal.setitimer(signal.ITIMER_REAL, timeout)

            with contextlib.redirect_stdout(self.stdout), contextlib.redirect_stderr(
                self.stderr
            ):
                try:
                    self.simulation.run(
                        steps=max(self.simulation.tsteps - self.simulation.step, 0)
                    )
                    self.returncode = 0
                except sp.TimeoutExpired:
                    raise
                except SystemExit as exc:
                    self.returncode = exc.code or 0
                except BaseException:
                    self.returncode = 1
                    self.stderr.write(traceback.format_exc())
                finally:
                    self._ran = True
                    self.kill()
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, old_alarm)
            os.chdir(old_cwd)

        return (
            self.stdout.getvalue().encode("ascii", errors="replace"),
            self.stderr.getvalue().encode("ascii", errors="replace"),
        )

    def kill(self):
        """Runs i-PI cleanup and resets soft-exit state.

        This mirrors the part of :meth:`subprocess.Popen.kill` that the
        regression-test runner needs: after it returns, sockets and callbacks
        are cleaned up and the wrapper has a non-zero return code unless the
        simulation had already set one.
        """

        if self.returncode is None:
            self.returncode = -9
        if self._cleanup_done:
            return
        old_cwd = os.getcwd()
        try:
            os.chdir(self.cwd)
            # Run softexit cleanup callbacks without calling sys.exit().
            with contextlib.redirect_stdout(self.stdout), contextlib.redirect_stderr(
                self.stderr
            ):
                softexit.cleanup()
        finally:
            os.chdir(old_cwd)
            self.simulation = None

        self._cleanup_done = True
        softexit.reset()


def get_test_settings(
    example_folder,
    settings_file="test_settings.dat",
):
    """This function looks for the existence of test_settings.dat file.
    This file can contain instructions like number of steps or driver name.
    If the file doesn't  exist, the driver dummy is assigned."""

    driver_models = list()
    driver_codes = list()
    socket_modes = list()
    port_numbers = list()
    address_names = list()
    flaglists = list()
    found_nsteps = False

    try:
        with open(Path(example_folder) / settings_file) as f:
            lines = f.readlines()
            if len(lines) == 0:
                raise ValueError("Error: The test_settings.dat file is empty.")

        starts = [i for i, line in enumerate(lines) if "driver_model" in line]

        for client in range(len(starts)):
            if client < len(starts) - 1:
                block = lines[starts[client] : starts[client + 1]]
            else:
                block = lines[starts[client] :]
            driver_code = "fortran"
            driver_model = "dummy"
            address_name = "localhost"
            port_number = 33333
            socket_mode = "unix"
            flaglist = {}
            for line in block:
                if "driver_code" in line:
                    driver_code = line.split()[1]
                elif "driver_model" in line:
                    driver_model = line.split()[1]
                elif "address" in line:
                    address_name = line.split()[1]
                elif "port" in line:
                    port_number = line.split()[1]
                elif "socket_mode" in line:
                    socket_mode = line.split()[1]
                elif "flags" in line:
                    flaglist = {line.split()[1]: line.split()[2:]}
                elif "nsteps" in line:
                    nsteps = line.split()[1]
                    found_nsteps = True

            # Checking that each driver has appropriate settings, if not, use default.
            if driver_code == "fortran" and driver_model not in fortran_driver_models:
                print(
                    "{} is not in fortran_driver_model list. We default to dummy".format(
                        driver_model
                    )
                )
                driver_model = "dummy"
            elif driver_code == "python" and driver_model not in python_driver_models:
                print(
                    "{} is not in python_driver_model list. We default to dummy".format(
                        driver_model
                    )
                )
                driver_model = "dummy"
            elif driver_code not in ["fortran", "python"]:
                raise ValueError(
                    "Driver code not available. Valid options are 'fortran' and 'python' only."
                )
            driver_models.append(driver_model)
            driver_codes.append(driver_code)
            socket_modes.append(socket_mode)
            port_numbers.append(port_number)
            address_names.append(address_name)
            flaglists.append(flaglist)

    except:
        driver_codes.append("fortran")
        driver_models.append("dummy")
        address_names.append("localhost")
        port_numbers.append(33333)
        socket_modes.append("unix")
        flaglists.append({})

    driver_info = {
        "driver_model": driver_models,
        "socket_mode": socket_modes,
        "port_number": port_numbers,
        "address_name": address_names,
        "driver_code": driver_codes,
        "flag": flaglists,
    }

    if found_nsteps:
        test_settings = {"nsteps": nsteps}
    else:
        test_settings = {"nsteps": "1"}

    return driver_info, test_settings


def modify_xml_4_dummy_test(
    input_name,
    output_name,
    nid,
    driver_info,
    test_settings,
):
    """Modify xml to run dummy tests"""
    try:
        tree = ET.parse(input_name)
    except:
        print("The error is in the format or the tags of the xml!")
    root = tree.getroot()
    clients = list()

    ff_roots = [root]
    if len(root.findall("ffcommittee")) > 0:
        ff_roots += root.findall("ffcommittee")
    if len(root.findall("ffrotations")) > 0:
        ff_roots += root.findall("ffrotations")

    ff_sockets = []
    for ff_root in ff_roots:
        for ff_socket in ff_root.findall("ffsocket"):
            ff_sockets.append(ff_socket)
        for ff_socket in ff_root.findall("ffcavphsocket"):
            ff_sockets.append(ff_socket)

    for s, ffsocket in enumerate(ff_sockets):
        ffsocket.attrib["mode"] = driver_info["socket_mode"][s]

        for element in ffsocket:
            port = driver_info["port_number"][s]
            if element.tag == "port":
                element.text = str(port)
            elif element.tag == "address":
                dd = driver_info["address_name"][s] + "_" + str(nid) + "_" + str(s)
                element.text = dd
                address = dd

        model = driver_info["driver_model"][s]

        clients.append([model, "unix", address, port])

        for key in driver_info["flag"][s].keys():
            if "-o" in key:
                for k, v in driver_info["flag"][s].items():
                    clients[s].append(k)
                    clients[s].extend(v)

        # if there are remaining drivers, assign them to the last socket
        if (s == len(ff_sockets) - 1) and (
            len(driver_info["driver_model"]) > len(ff_sockets)
        ):
            for remaining_client_idx in range(s + 1, len(driver_info["driver_model"])):
                model = driver_info["driver_model"][remaining_client_idx]
                clients.append([model, "unix", address, port])

                for key in driver_info["flag"][remaining_client_idx].keys():
                    if "-o" in key:
                        for k, v in driver_info["flag"][remaining_client_idx].items():
                            clients[remaining_client_idx].append(k)
                            clients[remaining_client_idx].extend(v)

    element = root.find("total_steps")
    if test_settings is not None:
        if element is not None:
            element.text = test_settings["nsteps"]
        else:
            new_element = ET.SubElement(root, "total_steps")
            new_element.text = test_settings["nsteps"]

    tree.write(open(output_name, "wb"))
    return clients


class Runner(object):
    """This class handles the creation of tmp directories,
    the i-pi call, the driver call, and finally
    it checks that the generated outputs are the expected ones.
    """

    def __init__(
        self,
        parent,
        call_ipi="i-pi input.xml",
        run_ipi_in_process=False,
    ):
        """Store parent directory and commands to call i-pi and driver
        call_ipi: command to call i-pi
        """

        self.parent = parent
        self.call_ipi = call_ipi
        self.run_ipi_in_process = run_ipi_in_process

    def run(self, cwd, nid):
        """This function tries to run the test in a tmp folder and
        checks if ipi has ended without error.
        Arguments:
            cwd: folder where all the original regression tests are stored
            nid: identification number to avoid repetitions of addresses"""

        try:
            # Create temp file and copy files
            self.tmp_dir = Path(tempfile.mkdtemp())
            print("\ntest folder: {}".format(cwd))
            print("temp folder: {}".format(self.tmp_dir))

            # files = os.listdir(self.parent / cwd)
            copy_tree(str(cwd), str(self.tmp_dir))
        except:
            return "Couldn't create the tmp folder"

        try:
            driver_info, test_settings = get_test_settings(self.tmp_dir)
            if "fortran" in driver_info["driver_code"]:
                call_driver = "i-pi-driver"
                address_key = "-h"
            elif "python" in driver_info["driver_code"]:
                call_driver = "i-pi-py_driver"
                address_key = "-a"
        except:
            print("Problem getting driver_info")
            return "Problem getting driver_info"

        clients = self.create_client_list(driver_info, nid, test_settings)

        try:
            # Run i-pi

            if getattr(self, "run_ipi_in_process", False):
                ipi = InProcessIPI(self.call_ipi, self.tmp_dir)
            else:
                ipi = sp.Popen(
                    self.call_ipi,
                    cwd=(self.tmp_dir),
                    shell=True,
                    stdin=sp.DEVNULL,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                )

            if len(clients) > 0:
                f_connected = False
                for client in clients:
                    for i in range(100):
                        if os.path.exists("/tmp/ipi_" + client[2]):
                            f_connected = True
                            break
                        else:
                            time.sleep(0.2)
                    if not f_connected:
                        # Check if i-pi finished successfully
                        if ipi.poll() is not None:
                            ipi_out, ipi_error = ipi.communicate(timeout=TIMEOUT)
                            if ipi.returncode == 0:
                                print(
                                    "i-PI finished before socket check. Assuming success."
                                )
                                return None
                            else:
                                print(
                                    "i-PI failed with return code {}".format(
                                        ipi.returncode
                                    )
                                )
                                print("i-PI Output:", ipi_out.decode("ascii"))
                                print("i-PI Error:", ipi_error.decode("ascii"))
                                return "i-PI failed with return code {}".format(
                                    ipi.returncode
                                )

                        print("Could not find the i-PI UNIX socket.")
                        print("Current client {}".format(client))
                        print("List all files  /tmp/ipi_*")
                        for filename in glob.glob("/tmp/ipi_*"):
                            print(filename)
                        ipi_out, ipi_error = ipi.communicate(timeout=TIMEOUT)
                        print("i-PI Output:", ipi_out.decode("ascii"))
                        print("i-PI Error:", ipi_error.decode("ascii"))
                        return "Could not find the i-PI UNIX socket"

            # Run drivers by defining cmd2 which will be called, eventually
            drivers = list()

            for client in clients:
                if client[1] == "unix":
                    clientcall = call_driver + " -m {} {} {} -u ".format(
                        client[0], address_key, client[2]
                    )
                elif client[1] == "inet":
                    clientcall = call_driver + " -m {} {} {} -p {}".format(
                        client[0], client[2], address_key, client[3]
                    )

                else:
                    raise ValueError("Driver mode has to be either unix or inet")

                cmd = clientcall

                # Add extra flags if necessary
                if any("-o" in str(s) for s in client):
                    flag_indeces = [
                        i for i, elem in enumerate(client) if "-o" in str(elem)
                    ]
                    for i, ll in enumerate(flag_indeces):
                        if i < len(flag_indeces) - 1:
                            cmd += " {} {}".format(
                                client[ll],
                                ",".join(client[ll + 1 : flag_indeces[i + 1]][:]),
                            )
                        else:
                            cmd += " {} {}".format(
                                client[ll], ",".join(client[ll + 1 :][:])
                            )

                if debug:
                    print(" Client call command:", cmd)

                driver = sp.Popen(
                    cmd,
                    cwd=(cwd),
                    shell=True,
                    stdin=sp.DEVNULL,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                )

                drivers.append(driver)

            # check i-pi errors
            ipi_out, ipi_error = ipi.communicate(timeout=TIMEOUT)
            assert ipi.returncode == 0, "i-PI error occurred: {}".format(ipi_error)

            # check driver errors
            for driver in drivers:
                # if i-PI has ended, we can wait for the driver to quit
                driver_out, driver_err = driver.communicate(timeout=TIMEOUT)
                if driver.returncode != 0:
                    ipi.kill()
                    ipi_out, ipi_error = ipi.communicate()
                    assert (
                        driver.returncode == 0
                    ), "Driver error occurred: {}\n Driver Output: {}\n i-PI Output: {}\n i-PI Error: {}".format(
                        driver_err, driver_out, ipi_out, ipi_error
                    )

        except sp.TimeoutExpired:
            ipi.kill()
            try:
                ipi_out, ipi_error = ipi.communicate(timeout=2)
            except:
                ipi_out, ipi_error = "", "Could not get outputs from ipi"
                pass

            drivers[0].kill()
            try:
                driver_out, driver_err = drivers[0].communicate(timeout=2)
            except:
                driver_out, driver_err = "", "Could not get outputs from drivers"
                pass

            print("Timeout during {} test \
              **** i-PI output **** \
              stdout {} \
              stderr {} \
              **** driver output **** \
              stdout {} \
              stderr {} \
              ".format(str(cwd), ipi_out, ipi_error, driver_out, driver_err))
            raise

            raise RuntimeError(
                "Time is out. Aborted during {} test.".format(
                    str(cwd),
                )
            )

        except FileNotFoundError:
            return "File not found\n{}".format(str(cwd))

        except ValueError:
            return "Value Error\n{}".format(str(cwd))
