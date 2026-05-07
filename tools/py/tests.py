#!/usr/bin/env python3

"""
This script runs calls pytest to run of the tests that are provided with i-pi.
For details of usage call it with "-h" option.
"""
import argparse
from argparse import RawTextHelpFormatter
from pathlib import Path
import pytest

import sys
import os

# Check that we have the import path for this i-PI set and if not, add it.
dir_root = os.path.realpath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "../..")
)
if dir_root not in sys.path:
    sys.path.insert(0, dir_root)

main_folder = Path(__file__).resolve().parents[2] / "ipi_tests"
test_folder = {
    "all": main_folder,
    "examples": main_folder / "examples",
    "regtests": main_folder / "regression_tests",
    "unit": main_folder / "unit_tests",
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description=""
        "Script that runs the tests provided with the code. \n"
        "\n"
        "If you want to run all the tests, \n"
        "type: i-pi-tests \n"
        "\n"
        "for running only the tests that checks the integrity of the examples inputs \n"
        "type: i-pi-tests -t examples  \n"
        "\n"
        "for running the regression tests  \n"
        "type: i-pi-tests -t regtests \n"
        "\n"
        "and for running the unitary tests  \n"
        "type: i-pi-tests -t unit \n",
    )

    parser.add_argument(
        "-t",
        "--tests",
        type=str,
        default="all",
        choices=["all", "examples", "regtests", "unit"],
        help="Specifies which tests are going to be called.",
    )

    args = parser.parse_args()
    tests = args.tests
    print(test_folder[args.tests])
    exit_code = pytest.main(
        ["--tb=line", "--capture=no", "-ra", "-v", str(test_folder[args.tests])]
    )

    if exit_code != 0:
        raise RuntimeError("pytest exit code is {}".format(exit_code))
