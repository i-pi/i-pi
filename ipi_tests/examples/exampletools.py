import subprocess as sp
import os
from pathlib import Path
from distutils.dir_util import copy_tree
import xml.etree.ElementTree as ET
import tempfile
import time

driver_models = [
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
    "doublewell",
    "doublewell_1D",
    "MB",
]


def get_single_driver_info(block, drivers, sockets, ports, addresses, flaglists):

    found_driver = False
    found_socket = False
    found_address = False
    found_port = False
    found_flags = False
    for line in block:
        if "driver" in line:
            driver = line.split()[1]
            if driver not in driver_models:
                driver = "dummy"
            drivers.append(driver)
            found_driver = True
        if "address" in line:
            addresses.append(line.split()[1])
            found_address = True
        if "port" in line:
            ports.append(line.split()[1])
            found_port = True
        if "socket_mode" in line:
            sockets.append(line.split()[1])
            found_socket = True
        if "flags" in line:
            flaglists.append({line.split()[1]: line.split()[2:]})
            found_flags = True

    # Checking that each driver has appropriate settings, if not, use default.
    if not found_driver:
        drivers.append("dummy")
    if not found_socket:
        sockets.append("unix")
    if not found_port:
        ports.append(33333)
    if not found_address:
        addresses.append("localhost")
    if not found_flags:
        flaglists.append({})

    return drivers, sockets, ports, addresses, flaglists


def get_test_settings(
    example_folder,
    settings_file="test_settings.dat",
    nsteps="2",
):
    """This function looks for the existence of test_settings.dat file.
    This file can contain instructions like number of steps or driver name.
    If the file doesn't  exist, the driver dummy is assigned."""
    drivers = list()
    sockets = list()
    ports = list()
    addresses = list()
    flaglists = list()

    try:
        with open(Path(example_folder) / settings_file) as f:
            lines = f.readlines()
            if len(lines) == 0:
                raise ValueError("test_settings.dat is empty")

        nline = 0
        starts = []
        for line in lines:
            if "driver" in line:
                starts.append(nline)
            nline += 1
            if "nsteps" in line:
                nsteps = line.split()[1]

        for client in range(len(starts)):
            if client < len(starts) - 1:
                block = lines[starts[client] : starts[client + 1]]
            else:
                block = lines[starts[client] :]
            drivers, sockets, ports, addresses, flaglists = get_single_driver_info(
                block, drivers, sockets, ports, addresses, flaglists
            )

    except:
        drivers, sockets, ports, addresses, flaglists = get_single_driver_info(
            "", drivers, sockets, ports, addresses, flaglists
        )
        pass

    driver_info = {
        "model": drivers,
        "socket_mode": sockets,
        "port_number": ports,
        "address_name": addresses,
        "flag": flaglists,
    }

    test_settings = {"nsteps": nsteps}

    return driver_info, test_settings


def find_examples(parent, excluded_file="excluded_test.txt", examples=[]):
    """This function looks for iteratively for examples and includes
    them if they don't appear in the excluded_file"""

    excluded = list()
    if excluded_file is not None:
        try:
            with open(excluded_file) as f:
                flines = [line for line in f.readlines() if line.strip()]
                for line in flines:
                    fname = "".join(str(parent / line).split()[0])
                    if line.split()[0] != "#":
                        excluded.append(fname)
        except:
            print("Excluded file not found")

    folders = [x[0] for x in os.walk(parent)]

    for ff in folders:
        if os.path.isfile(Path(ff) / "input.xml"):
            if ff not in excluded:
                examples.append(ff)

    return examples


def modify_xml_2_dummy_test(
    input_name,
    output_name,
    nid,
    driver_info,
    test_settings,
):
    """ Modify xml to run dummy tests """
    try:
        tree = ET.parse(input_name)
    except:
        print("The error is in the format or the tags of the xml!")
    root = tree.getroot()
    clients = list()

    if len(root.findall("ffcommittee")) > 0:
        ff_roots = root.findall("ffcommittee")
    else:
        ff_roots = [root]

    ff_sockets = []
    for ff_root in ff_roots:
        for ff_socket in ff_root.findall("ffsocket"):
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

        model = driver_info["model"][s]
        print("driver:", model)
        clients.append([model, address, port])

        for k, v in driver_info["flag"][s].items():
            clients[s].append(k)
            clients[s].extend(v)

    element = root.find("total_steps")
    if element is not None:
        element.text = test_settings["nsteps"]
    else:
        new_element = ET.SubElement(root, "total_steps")
        new_element.text = test_settings["nsteps"]

    tree.write(open(output_name, "wb"))
    return clients


class Runner_examples(object):
    """This class handles the modification of the examples inputs,
    the creation of tmp directories, the i-pi call, the driver call, and finally
    it checks if i-pi ended without error.
    """

    def __init__(self, parent, cmd1="i-pi new.xml", cmd2="i-pi-driver"):
        """ Store parent directory and commands to call i-pi and driver """
        self.parent = parent
        self.cmd1 = cmd1
        self.cmd2 = cmd2

    def run(self, cwd, nid):
        """This function tries to run the example in a tmp folder and
        afterwards checks if ipi has ended without error.
        arguments:
            cwd: folder where all the original examples files are stored
            nid: identification number to avoid repetitions of addresses"""

        try:
            # Create temp file and copy files
            self.tmp_dir = Path(tempfile.mkdtemp())
            print("temp folder: {}".format(self.tmp_dir))

            copy_tree(str(cwd), str(self.tmp_dir))
            driver_info, test_settings = get_test_settings(self.tmp_dir)

            # Modify xml
            clients = modify_xml_2_dummy_test(
                self.tmp_dir / "input.xml",
                self.tmp_dir / "new.xml",
                nid,
                driver_info,
                test_settings,
            )
            # Run i-pi
            ipi = sp.Popen(
                self.cmd1,
                cwd=(self.tmp_dir),
                shell=True,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
            )

            if len(clients) > 0:
                f_connected = False
                for i in range(50):
                    if os.path.exists("/tmp/ipi_" + clients[0][1]):
                        f_connected = True
                        break
                    else:
                        time.sleep(0.1)
                if not f_connected:
                    raise RuntimeError("Couldn't find the i-PI UNIX socket")

            # Run drivers
            driver = list()
            flag_indeces = list()
            for client in clients:
                cmd = self.cmd2 + " -m {} -u -h {}".format(client[0], client[1])
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

        except sp.TimeoutExpired:
            raise RuntimeError(
                "Time is out. Aborted during {} test. \
              Error {}".format(
                    str(cwd), ipi.communicate(timeout=2)[0]
                )
            )

        except FileNotFoundError:
            raise ("{}".format(str(cwd)))

        except ValueError:
            raise ("{}".format(str(cwd)))
