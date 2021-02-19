import subprocess as sp
import os
from pathlib import Path
from distutils.dir_util import copy_tree
import xml.etree.ElementTree as ET
import tempfile
import time
from ipi_tests.test_tools import get_test_settings, Runner


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


def modify_xml_4_dummy_test(
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

        clients.append([model, "unix", address, port])

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


class Runner_examples(Runner):
    """This class handles the modification of the examples inputs,
    the creation of tmp directories, the i-pi call, the driver call, and finally
    it checks if i-pi ended without error.
    """

    def __init__(self, parent, call_ipi="i-pi new.xml"):
        """ Store parent directory and commands to call i-pi """
        self.parent = parent
        self.call_ipi = call_ipi

    def create_client_list(self, driver_info, nid, test_settings):
        try:
            # Modify xml
            clients = modify_xml_4_dummy_test(
                self.tmp_dir / "input.xml",
                self.tmp_dir / "new.xml",
                nid,
                driver_info,
                test_settings,
            )
            return clients
        except:
            raise RuntimeError("Couldn't modify the xml file")
