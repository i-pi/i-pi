import subprocess as sp
from pathlib import Path
import pytest
import argparse
from argparse import RawTextHelpFormatter
import sys
from ipi_tests.examples.exampletools import find_examples, Runner_examples


""" Test that examples are not broken. Doesn't not check that output is correct."""

examples_folder = Path(__file__).resolve().parents[2] / "examples"
excluded_file = Path(__file__).resolve().parent / "excluded_test.txt"

examples = find_examples(examples_folder, excluded_file)
print("We have found {} reg_tests".format(len(examples)))


@pytest.mark.parametrize("ex", examples)
def test_example(ex):
    nid = examples.index(ex)
    runner = Runner_examples(Path("."))
    runner._run(ex, nid)


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
        "inclueded in the 'excluded_test.txt' file\n"
        "type: python test_examples.py \n"
        "\n"
        "To check all the examples inside a folder  \n"
        "type: python test_examples.py <folder_path> \n"
        "example: python test_examples.py examples/lammps/h2o-geop\n"
        "This script will recursively search for examples.\n",
    )

    parser.add_argument("folder", type=str, help="Folder of the example to test")
    args = parser.parse_args()
    try:
        path = Path(__file__).resolve().parents[2] / args.folder
        examples = find_examples(path, excluded_file=None)
        print("We will run only:")
        for i in examples:
            print(i)
    except:
        print("We will run all the tests")

    for ex in examples:
        print("Running {} ".format(ex))
        test_example(ex)
