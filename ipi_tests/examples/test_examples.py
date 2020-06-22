import subprocess as sp
from pathlib import Path
import pytest
import argparse
from argparse import RawTextHelpFormatter
import sys
import time
from ipi_tests.examples.exampletools import find_examples, Runner_examples


""" Test that examples are not broken. Doesn't not check that output is correct."""

examples_folder = Path(__file__).resolve().parents[2] / "examples"
excluded_file = Path(__file__).resolve().parent / "excluded_test.txt"

examples = find_examples(examples_folder, excluded_file)
print("We have found {} reg_tests".format(len(examples)))


@pytest.mark.parametrize("ex", examples)
def test_example(ex):
    t0 = time.time()
    nid = examples.index(ex)
    runner = Runner_examples(Path("."))
    runner.run(ex, nid)
    print("Time for this example: {:4.1f} s \n".format(time.time() - t0))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description=""
        "Script that check examples provided in the example folder\n"
        "It tests that i-pi runs and exists without error after few steps\n"
        "Doesn't not check that output is correct.\n"
        "\n"
        "It can be called by pytest or as a normal script \n"
        "To check all examples in the repository, except the ones \n"
        "included in the 'excluded_test.txt' file\n"
        "type: python test_examples.py \n"
        "\n"
        "To check all the examples inside a folder  except the ones \n"
        "inclueded in the 'excluded_test.txt' file\n"
        "type: python test_examples.py -f <folder_path> \n"
        "example: python test_examples.py examples/lammps/h2o-geop\n"
        "This script will recursively search for examples.\n"
        "\n"
        "To ignore excluded files use -test_all option",
    )

    parser.add_argument(
        "-f", "--folder", type=str, help="Folder of the example to test"
    )
    parser.add_argument(
        "--test_all", action="store_true", help="Folder of the example to test"
    )
    args = parser.parse_args()

    if args.test_all:
        excluded_file = None
    else:
        excluded_file = Path(__file__).resolve().parent / "excluded_test.txt"

    try:
        path = Path(__file__).resolve().parents[2] / args.folder
        examples = find_examples(path, excluded_file)
        print("We will run only:")
        for i in examples:
            print(i)
        print("")
    except:
        print("We will run all tests\n")
        examples = find_examples(examples_folder, excluded_file)

    errors = list()
    for ex in examples:
        print("Running {} ".format(ex))
        try:
            test_example(ex)
        except:
            print("Error ...\n")
            errors.append(ex)

    if len(errors) > 0:
        print("The following examples have exist with errors")
        for er in errors:
            print(er)
        print("")

    else:
        print("All the examples tested work.\n")
