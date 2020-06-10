import subprocess as sp
from pathlib import Path
import pytest
import time
import os
import shutil
import tempfile
from ipi_tests.regression.runstools import Runner, get_info_test
import argparse
from argparse import RawTextHelpFormatter


""" Run regression test """

main_folder = Path(__file__).parent
call_ipi = "i-pi input.xml"
call_driver = "i-pi-driver"

reg_tests = get_info_test(main_folder)


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

    parser = argparse.ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description=""
        "Script that performs regression tests\n"
        "It can be called by pytest or as a normal script \n"
        "To run all the regtest in the repository \n"
        "type: python test_run.py \n"
        "\n"
        "To check all the regtest inside a folder\n"
        "type: python test_run.py <folder_path> \n"
        "example: python test_run geop \n"
        "This script will recursively search for examples.\n",
    )

    parser.add_argument("folder", type=str, help="Folder of the example to test")
    args = parser.parse_args()

    path = main_folder / args.folder
    reg_tests = get_info_test(path)
    #try:
    #    path = main_folder / args.folder
    #    reg_tests = get_info_test(main_folder)
    #except:
    #    print("We will run all the tests")
    #    reg_tests = get_info_test(main_folder)


    print("We have found {} reg_tests".format(len(reg_tests)))
    for test_info in reg_tests:
        print("Running {} ".format(test_info[0]))
        test_cmd(test_info)
