import subprocess as sp
import os
from pathlib import Path
from distutils.dir_util import copy_tree
import xml.etree.ElementTree as ET
import tempfile
import time
from ipi_tests.test_tools import get_test_settings


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

    for s, ffsocket in enumerate(root.findall("ffsocket")):
        # name = ffsocket.attrib["name"]
        ffsocket.attrib["mode"] = driver_info["socket_mode"]

        for element in ffsocket:
            port = driver_info["port_number"]
            if element.tag == "port":
                element.text = str(port)
            elif element.tag == "address":
                dd = driver_info["address_name"] + "_" + str(nid) + "_" + str(s)
                element.text = dd
                address = dd

        model = driver_info["driver_model"]
        print("driver:", model)
        clients.append([model, address, port])

        for flag in driver_info["flag"]:
            for k, v in flag.items():
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

    def __init__(self, parent, cmd1="i-pi new.xml"):
        """ Store parent directory and commands to call i-pi """
        self.parent = parent
        self.cmd1 = cmd1

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
        except:
            return "Couldn't create the tmp folder"

        try:
            driver_info, test_settings = get_test_settings(self.tmp_dir)
            if driver_info["driver_code"] == "fortran":
                cmd2 = "i-pi-driver"
                address_key = "-h"
            elif driver_info["driver_code"] == "python":
                cmd2 = "i-pi-py_driver"
                address_key = "-a"
        except:
            return "Problem getting driver_info"

        try:
            # Modify xml
            clients = modify_xml_2_dummy_test(
                self.tmp_dir / "input.xml",
                self.tmp_dir / "new.xml",
                nid,
                driver_info,
                test_settings,
            )
        except:
            raise RuntimeError("Couldn't modify the xml file")
        print("Launching i-pi\n")

        try:
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
                    return "Couldn't find the i-PI UNIX socket"

            # Run drivers
            driver = list()
            flag_indeces = list()
            for client in clients:
                cmd = cmd2 + " -m {} -u {} {}".format(client[0], address_key, client[1])
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
