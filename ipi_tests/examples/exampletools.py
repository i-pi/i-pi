import os
from pathlib import Path

#!try:
#  from ..test_tools import Runner, modify_xml_4_dummy_test
# except:
from ipi_tests.test_tools import Runner, modify_xml_4_dummy_test


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
            if ff not in excluded and "broken" not in ff:
                examples.append(ff)

    return examples


class Runner_examples(Runner):
    """This class handles the modification of the examples inputs,
    the creation of tmp directories, the i-pi call, the driver call, and finally
    it checks if i-pi ended without error.
    """

    def __init__(self, parent, call_ipi="i-pi new.xml"):
        """Store parent directory and commands to call i-pi"""
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
        except Exception as ex:
            print("Couldn't modify the xml file")
            raise ex
