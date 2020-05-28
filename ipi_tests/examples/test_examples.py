import subprocess as sp
from pathlib import Path
import pytest
import time
import os
import sys
import shutil
import tempfile
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
    runner._run(ex,nid)


if __name__ == "__main__":
    try:
        examples = [sys.argv[1]]
        print("We will run only {} ".format(examples))
    except:
        print("We will run all the tests")

    for ex in examples:
        print("Running {} ".format(ex))
        test_example(ex)
