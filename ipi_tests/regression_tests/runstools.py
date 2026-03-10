import numpy as np
from pathlib import Path

try:
    from ..test_tools import (
        Runner,
        modify_xml_4_dummy_test,
    )
except ImportError:
    from ipi_tests.test_tools import Runner, modify_xml_4_dummy_test


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
        self.usetol = []

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
                    listcol = []
                    listtol = []
                    info_tol = line.split()[2].split(",")
                    if info_tol[0] == "all":
                        listcol = info_tol[0]
                        try:
                            listtol.append((float(info_tol[1]), float(info_tol[2])))
                        except:
                            listtol.append(None)
                    else:
                        # We have tol values for specific columns
                        for ll in line.split()[2:]:
                            lltol = ll.split(",")
                            listcol.append(int(lltol[0]))
                            if len(lltol) > 1:
                                listtol.append((float(lltol[1]), float(lltol[2])))
                            else:
                                listtol.append(None)
                    self.usecol.append(listcol)
                    self.usetol.append(listtol)
                else:
                    self.usecol.append(None)
                    self.usetol.append(None)

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
                    if type(self.usecol[ii]) is list:
                        usecol = self.usecol[ii]
                    else:
                        usecol = None
                    ref_output = np.loadtxt(Path(cwd) / refname, usecols=usecol)
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
                test_output = np.loadtxt(self.tmp_dir / fname, usecols=usecol)

                try:
                    if self.usetol[ii] is None:
                        rtol, atol = (1e-7, 1e-15)
                        np.testing.assert_allclose(
                            test_output, ref_output, rtol=rtol, atol=atol
                        )
                    else:
                        for icol, tol in enumerate(self.usetol[ii]):
                            if tol is None:
                                rtol, atol = (1e-7, 1e-15)
                            else:
                                rtol, atol = tol
                            np.testing.assert_allclose(
                                test_output[..., icol],
                                ref_output[..., icol],
                                rtol=rtol,
                                atol=atol,
                            )
                except AssertionError:
                    raise AssertionError(
                        "ANOMALY: Disagreement between reference and {} in {}. Absolute discrepancy (mean,max): {},{}.".format(
                            fname,
                            str((self.parent / cwd).absolute()),
                            np.abs(test_output - ref_output).mean(),
                            np.abs(test_output - ref_output).max(),
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
                        "Please provide a reference file named {} in {}\n \
                         (Note that extension *out  appears in gitignore\n \
                         so those files require to force the addition by\n \
                         'git add -f <filename.out>') ".format(
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
                    if self.usetol[ii] is None:
                        rtol, atol = (1e-7, 1e-6)
                        np.testing.assert_allclose(
                            test_xyz, ref_xyz, rtol=rtol, atol=atol
                        )
                    else:
                        for icol, tol in enumerate(self.usetol[ii]):
                            if tol is None:
                                rtol, atol = (1e-7, 1e-15)
                            else:
                                rtol, atol = tol
                            np.testing.assert_allclose(
                                test_xyz[..., icol],
                                ref_xyz[..., icol],
                                rtol=rtol,
                                atol=atol,
                            )
                except AssertionError:
                    raise AssertionError(
                        "ANOMALY: {} in {}. Absolute discrepancy (mean,max): {},{}.".format(
                            fname,
                            str((self.parent / cwd).absolute()).split("ipi_tests/", 1)[
                                1
                            ],
                            np.abs(test_xyz - ref_xyz).mean(),
                            np.abs(test_xyz - ref_xyz).max(),
                        )
                    )
