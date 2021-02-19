import subprocess as sp
import os
from pathlib import Path
from distutils.dir_util import copy_tree
import xml.etree.ElementTree as ET
import tempfile
import time

clean_all = False

fortran_driver_models = [
    "dummy",
    "lj",
    "sg",
    "harm",
    "harm3d",
    "morse",
    "zundel",
    "qtip4pf",
    "pswater",
    "eckart",
    "ch4hcbe",
    "MB",
]

# YL should do this  automatically but fFor now I do it explicitly
python_driver_models = ["dummy", "harmonic"]


def clean_tmp_dir():
    if clean_all:
        try:
            for filename in glob.glob("/tmp/ipi_*"):
                os.remove(filename)
        except:
            pass


def get_test_settings(
    example_folder,
    settings_file="test_settings.dat",
    driver_code="fortran",
    driver_model="dummy",
    socket_mode="unix",
    port_number=33333,
    address_name="localhost",
    flags=[],
    nsteps="2",
):
    """This function looks for the existence of test_settings.txt file.
    This file can contain instructions like number of steps or driver name.
    If the file doesn't  exist, the driver dummy is assigned."""
    try:
        with open(Path(example_folder) / settings_file) as f:
            flags = list()

            while True:
                line = f.readline()
                if not line:
                    break
                elif "driver_model" in line:
                    driver_model = line.split()[1]
                elif "socket_mode" in line:
                    socket_mode = line.split()[1]
                elif "port" in line:
                    port_number = line.split()[1]
                elif "address" in line:
                    address_name = line.split()[1]
                elif "flags" in line:
                    flags.append({line.split()[1]: line.split()[2:]})
                elif "nsteps" in line:
                    nsteps = line.split()[1]
                elif "driver_code" in line:
                    driver_code = line.split()[1]
    except:
        pass

    if driver_code == "fortran":
        if driver_model not in fortran_driver_models:
            driver_model = "dummy"
    elif driver_code == "python":
        if driver_model not in python_driver_models:
            driver_model = "dummy"
    else:
        raise ValueError(
            "Drive code not available. Valid options are 'fortran' and 'python'"
        )

    driver_info = {
        "driver_model": driver_model,
        "socket_mode": socket_mode,
        "port_number": port_number,
        "address_name": address_name,
        "driver_code": driver_code,
        "flag": flags,
    }

    test_settings = {"nsteps": nsteps}

    return driver_info, test_settings


class Runner(object):
    """This class handles the creation of tmp directories,
    the i-pi call, the driver call, and finally
    it checks that the generated outputs are the expected ones.
    """

    def __init__(
        self,
        parent,
        call_ipi="i-pi input.xml",
    ):
        """Store parent directory and commands to call i-pi and driver
        call_ipi: command to call i-pi
        """

        self.parent = parent
        self.call_ipi = call_ipi

    def run(self, cwd, nid):
        """This function tries to run the test in a tmp folder and
        checks if ipi has ended without error.
        Arguments:
            cwd: folder where all the original regression tests are stored
            nid: identification number to avoid repetitions of addresses"""

        try:
            # Create temp file and copy files
            self.tmp_dir = Path(tempfile.mkdtemp())
            print("\nTest folder: {}".format(cwd))
            print("temp folder: {}".format(self.tmp_dir))

            files = os.listdir(self.parent / cwd)
            copy_tree(str(cwd), str(self.tmp_dir))
        except:
            return "Couldn't create the tmp folder"

        try:
            driver_info, test_settings = get_test_settings(self.tmp_dir)
            if driver_info["driver_code"] == "fortran":
                call_driver = "i-pi-driver"
                address_key = "-h"
            elif driver_info["driver_code"] == "python":
                call_driver = "i-pi-py_driver"
                address_key = "-a"
        except:
            return "Problem getting driver_info"

        clients = self.create_client_list(driver_info, nid, test_settings)

        if True:
            # try:
            # Run i-pi

            ipi = sp.Popen(
                self.call_ipi,
                cwd=(self.tmp_dir),
                shell=True,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
            )

            if len(clients) > 0:
                f_connected = False
                for i in range(50):
                    if os.path.exists("/tmp/ipi_" + clients[0][2]):
                        f_connected = True
                        break
                    else:
                        time.sleep(0.1)
                if not f_connected:
                    return "Couldn't find the i-PI UNIX socket"

            # Run drivers by defining cmd2 which will be called, eventually
            driver = list()
            flag_indeces = list()

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
                if any("-" in str(s) for s in client):
                    flag_indeces = [
                        i for i, elem in enumerate(client) if "-" in str(elem)
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
                print("cmd:", cmd)
                driver.append(
                    sp.Popen(cmd, cwd=(cwd), shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
                )

            # Check errors
            ipi_error = ipi.communicate(timeout=120)[1].decode("ascii")
            if ipi_error != "":
                print(ipi_error)
            assert "" == ipi_error

        try:
            2 == 2
        except sp.TimeoutExpired:
            raise RuntimeError(
                "Time is out. Aborted during {} test. \
              Error {}".format(
                    str(cwd), ipi.communicate(timeout=2)[0]
                )
            )

        except FileNotFoundError:
            return "File not found\n{}".format(str(cwd))

        except ValueError:
            return "Value Error\n{}".format(str(cwd))
