import subprocess as sp
import os
from pathlib import Path
import xml.etree.ElementTree as ET
import sys
import tempfile
import shutil
import time

driver_models = [
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
]


def find_examples(parent, excluded_file="excluded_test.txt"):

    excluded = list()
    with open(excluded_file) as f:
        for line in f:
            excluded.append(" ".join(str(parent / line).split()))

    folders = [x[0] for x in os.walk(parent)]
    examples = list()

    for ff in folders:
        if os.path.isfile(Path(ff) / "input.xml"):
            if ff not in excluded:
                examples.append(ff)
    return examples


def modify_xml(
    input_name,
    output_name,
    socket_mode="unix",
    port_number=33333,
    address_name="localhost",
    nsteps=5,
):
    """ Modify xml to run dummy tests """

    tree = ET.parse(input_name)
    root = tree.getroot()
    clients = list()

    for ffsocket in root.findall("ffsocket"):
        name = ffsocket.attrib["name"]
        ffsocket.attrib["mode"] = socket_mode

        for element in ffsocket:
            port = port_number + len(clients)
            if element.tag == "port":
                element.text = str(port)
            elif element.tag == "address":
                element.text = address_name + "_" + str(port)
                address = address_name + "_" + str(port)
        model = ffsocket.attrib["name"]
        if model not in driver_models:
            model = "gas"
        clients.append((model, address, port))

    element = root.find("total_steps")
    if element is not None:
        element.text = str(nsteps)
    else:
        new_element = ET.SubElement(root, "total_steps")
        new_element.text = str(nsteps)

    tree.write(open(output_name, "wb"))
    return clients


class Runner_examples(object):
    def __init__(self, parent, cmd1="i-pi new.xml", cmd2="i-pi-driver"):
        self.parent = parent
        self.cmd1 = cmd1
        self.cmd2 = cmd2

    def _run(self, cwd):
        try:
            self.tmp_dir = Path(tempfile.mkdtemp())
            files = os.listdir(self.parent / cwd)
            for f in files:
                shutil.copy(self.parent / cwd / f, self.tmp_dir)
            clients = modify_xml(
                self.tmp_dir / "input.xml", self.tmp_dir / "new.xml", socket_mode="unix"
            )

            ipi = sp.Popen(
                self.cmd1,
                cwd=(self.tmp_dir),
                shell=True,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
            )
            time.sleep(2)

            driver = list()
            for client in clients:
                cmd = self.cmd2 + " -m {} -h {} -u ".format(client[0], client[1])
                driver.append(sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE))

            ipi_error = ipi.communicate(timeout=30)[1].decode("ascii")
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


if __name__ == "__main__":
    name = sys.argv[1]
    out = "output.xml"
    modify_xml(name, out)
