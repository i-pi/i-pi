import subprocess as sp
import os
import numpy as np
from pathlib import Path
import xml.etree.ElementTree as ET
import tempfile
import shutil
import time
import glob

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
    "ljpolymer",
    "doublewell",
    "doublewell_1D",
    "MB",
]


def get_driver_info(
    example_folder,
    driver_info_file="driver.txt",
    driver="gas",
    socket_mode="unix",
    port_number=33333,
    address_name="localhost",
    flags=[],
):
    """This function looks for the existence of .driver_info file
    to run the example with a meaningfull driver. If the file doesn't
    exist, the driver gas is assigned."""
    try:
        with open(Path(example_folder) / driver_info_file) as f:
            flags = list()
            ncount = 0
            while True:
                line = f.readline()
                if not line:
                    break
                elif "driver" in line:
                    driver = line.split()[1]
                elif "socket_mode" in line:
                    socket_mode = line.split()[1]
                elif "port" in line:
                    port_number = line.split()[1]
                elif "address" in line:
                    address_name = line.split()[1]
                elif "flags" in line:
                    flags.append({line.split()[1]: line.split()[2:]})
                ncount += 1
            if ncount == 0:
                raise ValueError("driver.txt is empty")
    except FileNotFoundError:
        raise FileNotFoundError(
            "({}) An input.xml file was found but a driver.txt not".format(
                example_folder
            )
        )
    except:
        pass

    if driver not in driver_models:
        driver = "gas"

    driver_info = {
        "model": driver,
        "socket_mode": socket_mode,
        "port_number": port_number,
        "address_name": address_name,
        "flag": flags,
    }

    return driver_info


def get_info_test(parent):
    """This function recursively searches for examples
    and checks for the presence of the required additional info
    required (i.e. driver.txt)
    """
    folders = [x[0] for x in os.walk(parent)]
    reg_tests = list()

    for ff in folders:
        if os.path.isfile(Path(ff) / "input.xml"):
            try:
                driver_info = get_driver_info(ff)
                reg_tests.append([ff, driver_info])
            except FileNotFoundError:
                raise FileNotFoundError(
                    "({}) An input.xml file was found but a driver.txt not".format(ff)
                )
            except ValueError:
                raise ValueError(
                    'Please specify "model address      port mode" inside driver.txt'
                )

    return reg_tests


class Runner(object):
    """This class handles the creation of tmp directories,
    the i-pi call, the driver call, and finally
    it checks that the generated outputs are the expected ones.
    """

    def __init__(
        self,
        parent,
        call_ipi="i-pi input.xml",
        call_driver="i-pi-driver",
        check_errors=True,
        check_numpy_output=True,
        check_xyz_output=True,
    ):
        """Store parent directory and commands to call i-pi and driver
        call_ipi: command to call i-pi
        call_driver: list of commands to call drivers
        """

        self.parent = parent
        self.call_ipi = call_ipi
        self.call_driver = call_driver
        self.check_error = check_errors
        self.check_numpy_output = check_numpy_output
        self.check_xyz_output = check_xyz_output

    def _run(self, info, nid):
        """This function tries to run the example in a tmp folder and
        checks if ipi has ended without error.
        After that the output is checked against a reference
        arguments:
            cwd: folder where all the original regression tests are stored
        """

        cwd = info[0]
        self.files = []
        self.forms = []

        try:
            # Create temp file and copy files
            self.tmp_dir = Path(tempfile.mkdtemp())
            print("temp folder: {}".format(self.tmp_dir))

            files = os.listdir(self.parent / cwd)
            for f in files:
                shutil.copy(self.parent / cwd / f, self.tmp_dir)
            driver_info = get_driver_info(self.tmp_dir)

            # Creating the clients list
            input_name = self.tmp_dir / "input.xml"
            output_name = self.tmp_dir / "input.xml"
            try:
                tree = ET.parse(input_name)
            except:
                print("The error is in the format or the tags of the xml!")
            root = tree.getroot()

            clients = list()

            for s, ffsocket in enumerate(root.findall("ffsocket")):
                #name = ffsocket.attrib["name"]
                mode = driver_info["socket_mode"]
                ffsocket.attrib["mode"] = mode

                for element in ffsocket:
                    port = driver_info["port_number"]
                    if element.tag == "port":
                        element.text = str(port)
                    elif element.tag == "address":
                        dd = driver_info["address_name"] + "_" + str(nid) + "_" + str(s)
                        element.text = dd
                        address = dd

                model = driver_info["model"]
                clients.append([model, address, port, mode])

                for flag in driver_info["flag"]:
                    for k, v in flag.items():
                        clients[s].append(k)
                        clients[s].extend(v)

            tree.write(open(output_name, "wb"))

            # Run i-pi
            ipi = sp.Popen(
                self.call_ipi,
                cwd=(self.tmp_dir),
                shell=True,
                stdout=sp.PIPE,
                stderr=sp.PIPE,
            )
            time.sleep(2)

            with open(self.tmp_dir / "files_to_check.txt") as f:
                lines = f.readlines()
            for nl, line in enumerate(lines):
                if nl > 1:
                    self.files.append(line.split()[0])
                    self.forms.append(line.split()[1])

            # Run drivers by defining cmd2 which will be called, eventually
            flag_indeces = list()

            driver = list()
            cmd2 = list()
            for client in clients:
                if client[3] == "unix":
                    clientcall = self.call_driver + " -m {} -h {} -u ".format(
                        client[0], client[1]
                    )
                    cmd2.append(clientcall)
                elif client[3] == "inet":
                    cmd2.append(
                        self.call_driver
                        + " -m {} -h {} -p {}".format(client[0], client[1], client[2])
                    )
                else:
                    raise ValueError("Driver mode has to be either unix or inet")

                cmd = cmd2[0]
                if any("-" in str(s) for s in client):
                    flag_indeces = [
                        i for i, elem in enumerate(client) if "-" in str(elem)
                    ]
                    for i, ll in enumerate(flag_indeces):
                        if i < len(flag_indeces) - 1:
                            cmd2 += "{} {}".format(
                                client[ll],
                                ",".join(client[ll + 1 : flag_indeces[i + 1]][:]),
                            )
                        else:
                            cmd2.append(
                                " {} {}".format(
                                    client[ll], ",".join(client[ll + 1 :][:])
                                )
                            )
                    cmd += cmd2[1]
                driver.append(
                    sp.call(cmd, cwd=(cwd), shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
                )

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
        try:
            for filename in glob.glob("/tmp/ipi_*"):
                os.remove(filename)
        except:
            pass

        if self.check_error:
            self._check_error(ipi)
        if self.check_numpy_output:
            self._check_numpy_output(cwd)
        if self.check_xyz_output:
            self._check_xyz_output(cwd)

    def _check_error(self, ipi):
        """ This function checks if ipi has exited with errors"""

        ipi_error = ipi.communicate(timeout=120)[1].decode("ascii")
        if ipi_error != "":
            print("IPI ERROR OCCURED: {}".format(ipi_error))
        assert "" == ipi_error

    def _check_numpy_output(self, cwd):
        """This function checks if the numpy-accessible datafiles are 'all_close' to the
        reference file provided

        The checked 'numpy-accessible' files are collected in the files_to_check.txt
        where there filetype is specified, too.
        """

        for ii, refname in enumerate(self.files):
            if self.forms[ii] == "numpy":
                try:
                    ref_output = np.loadtxt(Path(cwd) / refname)
                except IOError:
                    raise IOError(
                        'Please provide a reference properties output named "{}"'.format(
                            refname
                        )
                    )
                except ValueError:
                    raise ValueError(
                        "Please check ref_simulation.out in {}".format(str(self.parent))
                    )

                fname = refname[4:]
                test_output = np.loadtxt(self.tmp_dir / fname)

                try:
                    np.testing.assert_allclose(
                        test_output, ref_output, rtol=1.0e-7, atol=1.0e-15
                    )
                    print("No anomaly during the regtest for {}".format(refname))
                except AssertionError:
                    raise AssertionError(
                        "ANOMALY: Disagreement between reference and {} in {}".format(
                            fname, str(self.parent)
                        )
                    )

    def _check_xyz_output(self, cwd):
        """This function checks if the ref_simulation.XXXXX.xyz files are 'all_close'
        to the reference file provided.

        The checked .xyz files are collected in the files_to_check.txt where there filetype is
        also specified.

        Common options for XXXXX and the tags they can be generated with:
            pos_c = <trajectory filename='pos_c' stride='1' > x_centroid </trajectory>
            frc_c = <trajectory filename='frc_c' stride='1' > f_centroid </trajectory>
            vel_c = <trajectory filename='vel_c' stride='1' > v_centroid </trajectory>
            mom_c = <trajectory filename='mom_c' stride='1' > p_centroid </trajectory>

        In case of a multi-bead simulations, please check of the positions and
        forces on the first bead, too.
            pos_0 = <trajectory filename='pos' stride='1' bead='0'> positions </trajectory>
            frc_0 = <trajectory filename='frc' stride='1' bead='0'> forces </trajectory>
            vel_0 = <trajectory filename='vel' stride='1' bead='0'> velocities </trajectory>
            mom_0 = <trajectory filename='mom' stride='1' bead='0'> momenta </trajectory>

        """

        for ii, refname in enumerate(self.files):
            if self.forms[ii] == "xyz":
                ref_structs = []
                try:
                    with open(Path(cwd) / refname) as ref:
                        ref_natoms = int(ref.readline())
                        for s, ll in enumerate(ref.readlines()):
                            if (s + 1) % (ref_natoms + 2) != 0 and (s + 1) % (
                                ref_natoms + 2
                            ) != 1:
                                ref_structs.append(ll.split()[1:])
                    reff = [[float(v) for v in r] for r in ref_structs]
                    ref_xyz = np.array(reff)
                except IOError:
                    raise IOError(
                        "Please provide a reference file named {} in {}".format(
                            refname, str(self.parent / cwd)
                        )
                    )

                except ValueError:
                    raise ValueError(
                        "Please check the values for the file named {} in {}".format(
                            refname, str(self.parent / cwd)
                        )
                    )

                fname = refname[4:]
                structs = []
                with open(self.tmp_dir / fname) as f:
                    natoms = int(f.readline())
                    for s, ll in enumerate(f.readlines()):
                        if (s + 1) % (natoms + 2) != 0 and (s + 1) % (natoms + 2) != 1:
                            structs.append(ll.split()[1:])
                    testt = [[float(v) for v in r] for r in structs]
                    test_xyz = np.array(testt)

                try:
                    np.testing.assert_allclose(
                        test_xyz, ref_xyz, rtol=1.0e-7, atol=1.0e-8
                    )
                    print("No anomaly during the regtest for {}".format(refname))
                except AssertionError:
                    raise AssertionError(
                        "ANOMALY: Disagreement between reference and {} in {}".format(
                            fname, str(self.parent / cwd)
                        )
                    )
