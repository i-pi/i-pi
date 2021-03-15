import numpy as np
from pathlib import Path
from ..test_tools import (
    Runner,
    modify_xml_4_dummy_test,
)


class Runner_regression(Runner):
    """This class handles the creation of tmp directories,
    the i-pi call, the driver call, and finally
    it checks that the generated outputs are the expected ones.
    """

    def __init__(
        self,
        parent,
        call_ipi="i-pi input.xml",
        check_numpy_output=True,
        check_xyz_output=True,
    ):
        """Store parent directory and commands to call i-pi and driver
        call_ipi: command to call i-pi
        call_driver: list of commands to call drivers
        """

        self.parent = parent
        self.call_ipi = call_ipi
        self.check_numpy_output = check_numpy_output
        self.check_xyz_output = check_xyz_output

        self.files = []
        self.forms = []
        self.usecol = []

    def create_client_list(self, driver_info, nid, test_settings):

        try:
            # Modify xml
            clients = modify_xml_4_dummy_test(
                self.tmp_dir / "input.xml",
                self.tmp_dir / "input.xml",
                nid,
                driver_info,
                test_settings=None,
            )
            return clients
        except:
            raise RuntimeError("Couldn't modify the xml file")

    def run(self, cwd, nid):
        """This function tries to run the example in a tmp folder and
        checks if ipi has ended without error.
        After that the output is checked against a reference
        arguments:
            cwd: folder where all the original regression tests are stored
        """

        super().run(cwd, nid)

        with open(self.tmp_dir / "files_to_check.txt") as f:
            lines = f.readlines()
        for nl, line in enumerate(lines):
            if nl > 1:
                self.files.append(line.split()[0])
                self.forms.append(line.split()[1])
                if len(line.split()) > 2:
                    listll = []
                    for ll in line.split()[2:]:
                        listll.append(int(ll))
                    self.usecol.append(listll)
                else:
                    self.usecol.append(None)

        if self.check_numpy_output:
            self._check_numpy_output(cwd)
        if self.check_xyz_output:
            self._check_xyz_output(cwd)

    def _check_numpy_output(self, cwd):
        """This function checks if the numpy-accessible datafiles are 'all_close' to the
        reference file provided

        The checked 'numpy-accessible' files are collected in the files_to_check.txt
        where there filetype is specified, too.
        """

        for ii, refname in enumerate(self.files):
            if self.forms[ii] == "numpy":
                try:
                    ref_output = np.loadtxt(
                        Path(cwd) / refname, usecols=self.usecol[ii]
                    )
                except IOError:
                    raise IOError(
                        'Please provide a reference properties output named "{}"'.format(
                            refname
                        )
                    )
                except ValueError:
                    raise ValueError(
                        "Please check ref_simulation.out in {}".format(
                            str((self.parent / cwd).absolute())
                        )
                    )

                fname = refname[4:]
                test_output = np.loadtxt(self.tmp_dir / fname, usecols=self.usecol[ii])

                try:
                    np.testing.assert_allclose(
                        test_output, ref_output, rtol=1.0e-7, atol=1.0e-15
                    )
                    # print("No anomaly during the regtest for {}".format(refname))
                except AssertionError:
                    raise AssertionError(
                        "ANOMALY: Disagreement between reference and {} in {}".format(
                            fname, str((self.parent / cwd).absolute())
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
                            refname, str((self.parent / cwd).absolute())
                        )
                    )

                except ValueError:
                    raise ValueError(
                        "Please check the values for the file named {} in {}".format(
                            refname, str((self.parent / cwd).absolute())
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
                    # print("No anomaly during the regtest for {}".format(refname))
                except AssertionError:
                    raise AssertionError(
                        "ANOMALY: {} in {}".format(
                            fname,
                            str((self.parent / cwd).absolute()).split("ipi_tests/", 1)[
                                1
                            ],
                        )
                    )
