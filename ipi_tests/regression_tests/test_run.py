import subprocess as sp
from pathlib import Path
import pytest
import argparse
from argparse import RawTextHelpFormatter
import time
from ipi_tests.regression_tests.runstools import Runner, get_info_test


""" Run regression test """

regtests_folder = Path(__file__).resolve().parent / "tests"
call_ipi = "i-pi input.xml"
call_driver = "i-pi-driver"

reg_tests = get_info_test(regtests_folder)


@pytest.mark.parametrize("regtest", reg_tests)
def test_regtest(regtest):
    """ Intermediate function to run the regression test (by calling Runner) and makes
    possible to parametrize the arguments
    """
    t0 = time.time()
    nid = reg_tests.index(regtest)
    # print('nid=',nid)
    runner = Runner(Path("."))

    cmd2 = list()
    # for t in regtest[1]:
    if regtest[1]["socket_mode"] == "unix":
        cmd2.append(
            call_driver
            + " -m {} -h {} -u ".format(regtest[1]["model"], regtest[1]["address_name"])
        )
    elif regtest[1]["socket_mode"] == "inet":
        cmd2.append(
            call_driver
            + " -m {} -h {} -p {}".format(
                regtest[1]["model"],
                regtest[1]["address_name"],
                regtest[1]["port_number"],
            )
        )
    else:
        raise ValueError("Driver mode has to be either unix or inet")

    runner._run(test_info[0], nid)
    print("Time for this regtest: {:4.1f} s \n".format(time.time() - t0))


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
        "type: python test_run.py --path <folder_path> \n"
        "example: python test_run -p geop \n"
        "This script will recursively search for examples.\n",
    )

    parser.add_argument(
        "-f",
        "--folder",
        type=str,
        default=None,
        help="Folder of the regressions to test",
    )
    parser.add_argument(
        "--test_all",
        action="store_true",
        help="Shall we test all of the regression-examples ?",
    )
    args = parser.parse_args()

    try:
        path = regtests_folder / args.folder
        reg_tests = get_info_test(path)
        print("We will run only:")
        for i in reg_tests:
            print(i[0])
        print("")
    except:
        print("We will run all available regression tests")
        reg_tests = get_info_test(regtests_folder)

    print("We have found {} reg_tests".format(len(reg_tests)))
    for test_info in reg_tests:
        print("Running {} ".format(test_info[0]))
        test_regtest(test_info)
