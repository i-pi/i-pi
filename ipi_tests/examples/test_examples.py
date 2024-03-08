from pathlib import Path
import pytest
import argparse
from argparse import RawTextHelpFormatter
import time

try:
    from .exampletools import find_examples, Runner_examples
except:
    from exampletools import find_examples, Runner_examples

""" Test that examples are not broken. Doesn't not check that output is correct."""

examples_folder = Path(__file__).resolve().parents[2] / "examples"
excluded_file = Path(__file__).resolve().parent / "excluded_test.txt"

examples = []
examples = find_examples(examples_folder, excluded_file, examples)
if __name__ != "__main__":
    print("We have found {} examples".format(len(examples)))


@pytest.mark.parametrize("ex", examples)
def test_example(ex, verbose=False):
    """Intermediate function to run the examples (by calling Runner_examples) which makes
    possible to parametrize the arguments
    """
    t0 = time.time()
    nid = examples.index(ex)
    runner = Runner_examples(Path("."))
    error_msg = runner.run(ex, nid)
    print("Time for this example: {:4.1f} s \n".format(time.time() - t0))

    if verbose:
        return error_msg

    if error_msg != None:
        raise RuntimeError


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
        "\n"
        "type: python test_examples.py \n"
        "\n"
        "\n"
        "To check all the examples inside a folder  except the ones \n"
        "included in the 'excluded_test.txt' file\n"
        "\n"
        "type: python test_examples.py -f <folder_path> \n"
        "example  python test_examples.py -f examples/features/ring_polymer_instanton examples/hpc_scripts/slurm  \n"
        "Note that the folder path is referenced to the i-pi root folder  \n"
        "\n"
        "example: python test_examples.py examples/lammps/h2o-geop\n"
        "This script will recursively search for examples.\n"
        "\n"
        "To ignore excluded files use -test_all option",
    )

    parser.add_argument(
        "-f",
        "--folder",
        type=str,
        nargs="+",
        help="Folder(s) of the example to test. Do not give absolute path. The path should be referenced to the i-pi root directory. Example: -f examples/MBPOL/splitting",
    )
    parser.add_argument(
        "--test_all", action="store_true", help="Folder of the example to test"
    )
    args = parser.parse_args()

    if args.test_all:
        excluded_file = None
    else:
        excluded_file = Path(__file__).resolve().parent / "excluded_test.txt"

    examples = list()
    try:
        for folder in args.folder:
            path = Path(__file__).resolve().parents[2] / folder
            print("Looking for examples inside {}".format(str(path)))
            examples.extend(find_examples(path, excluded_file, examples))
        print("\nWe will run tests in:")
        for i in examples:
            print(i)
        print("")
    except:
        print("We will run all tests\n")
        examples = find_examples(examples_folder, excluded_file, examples)

    errors = list()
    for ex in examples:
        print("Running {} ".format(ex))
        error_msg = test_example(ex, verbose=True)
        if error_msg is not None:
            print("Error: {}...\n".format(error_msg))
            errors.append(ex)

    if len(errors) > 0:
        print("The following examples have exited with errors")
        for er in errors:
            print(er)
        print("")

    else:
        print("All the examples tested work.\n")
