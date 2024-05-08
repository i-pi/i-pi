from pathlib import Path
import pytest
import argparse
from argparse import RawTextHelpFormatter
import time

try:
    from ipi_tests.regression_tests.runstools import Runner_regression
    from ipi_tests.test_tools import get_test_list
except ImportError:
    from runstools import Runner_regression
    from ipi_tests.test_tools import get_test_list


""" Run regression test """

regtests_folder = Path(__file__).resolve().parent / "tests"
call_ipi = "i-pi input.xml"

reg_tests = get_test_list(regtests_folder)


@pytest.mark.parametrize("regtest", reg_tests)
def test_regtest(regtest):
    """Intermediate function to run the regression test (by calling Runner) and makes
    possible to parametrize the arguments
    """
    t0 = time.time()
    nid = reg_tests.index(regtest)
    runner = Runner_regression(Path("."))
    error_msg = runner.run(regtest, nid)
    print("Time for this example: {:4.1f} s \n".format(time.time() - t0))

    if error_msg != None:
        raise RuntimeError(error_msg)


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
        help="Folder of the regressions to test. Example '-f GEOP'",
    )
    parser.add_argument(
        "--test_all",
        action="store_true",
        help="Shall we test all of the regression-examples ?",
    )
    args = parser.parse_args()

    try:
        path = str(regtests_folder / args.folder)
        reg_tests = get_test_list(path, skip_noauto=False)
        print("We will run only:")
        for i in reg_tests:
            print(i)
        print("")
    except:
        print("We will run all available regression tests")
        reg_tests = get_test_list(regtests_folder, skip_noauto=True)

    print("We have found {} reg_tests".format(len(reg_tests)))
    for test_info in reg_tests:
        print("Running {} ".format(test_info))
        test_regtest(test_info)
