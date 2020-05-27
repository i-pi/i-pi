import subprocess as sp
from pathlib import Path
import pytest
import time
import os
import shutil
import tempfile
from ipi_tests.reg_test.runstools import Runner, get_info_test

""" Run regression test """

main_folder = Path(__file__).parent
call_ipi = "i-pi input.xml"
call_driver = "i-pi-driver"

reg_tests = get_info_test(main_folder)
print("We have found {} reg_tests".format(len(reg_tests)))


@pytest.mark.parametrize("test_info", reg_tests)
def test_cmd(test_info):
    runner = Runner(Path("."))

    cmd2 = list()
    for t in test_info[1]:
        if t[3] == "unix":
            cmd2.append(call_driver + " -m {} -h {} -u ".format(t[0], t[1]))
        elif t[3] == "inet":
            cmd2.append(call_driver + " -m {} -h {} -p {}".format(t[0], t[1], t[2]))
        else:
            raise ValueError("Driver mode has to be either unix or inet")
    runner._run(call_ipi, cmd2, test_info[0])


if __name__ == "__main__":

    for test_info in reg_tests:
        print("Running {} ".format(test_info[0]))
        test_cmd(test_info)
