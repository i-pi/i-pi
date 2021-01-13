#!/usr/bin/env python3

"""
This script cleans python files in-place by calling "formatter" variable
and checks linting by using "linter" variable.
Prior to running the script, black >20.8b1  must be installed.
For details of usage call it with "-h" option.
"""
import subprocess as sp
import argparse
from argparse import RawTextHelpFormatter
from pathlib import Path
import sys

formatter = "black"
linter = "black --check"
cwd_repo = Path("./").resolve()


def main(check_mode, folder_path, file_path):

    if file_path is not None:
        path = cwd_repo / file_path
    elif folder_path is not None:
        path = cwd_repo / folder_path
    else:
        path = cwd_repo

    if not check_mode:
        print("\nTarget: {}".format(path))
        cmd1 = [formatter, path]
        print("\nRun {}".format(cmd1))
        formatt = sp.Popen(cmd1, cwd=cwd_repo, stdout=sp.PIPE, stderr=sp.PIPE)

        msg = formatt.communicate(timeout=30)[1].decode("utf-8")
        print("{} output:\n".format(formatter))
        print(msg)

    print("\nRun {}".format(linter))
    cmd1 = linter.split()
    cmd1.append(path)
    lint = sp.Popen(cmd1, cwd=cwd_repo, stdout=sp.PIPE, stderr=sp.PIPE)
    msg = lint.communicate(timeout=30)[0].decode("utf-8")
    if msg == "":
        print("\nThere isn't any warning. You are ready to push this commit.\n")
    else:
        print("{} output:\n".format(linter))
        print(msg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description=""
        "Script that takes care of the style guide enforcement. \n"
        "To run this script black >20.8b1 must be installed: \n"
        "https://black.readthedocs.io/en/stable/  \n"
        "\n"
        "If you want to clean your i-pi repository, \n"
        "type: i-pi-style \n"
        "\n"
        "If you want to clean a folder and everything inside it \n"
        "type: i-pi-style --path <folder_path> \n"
        "example: i-pi-style -p ipi/engine/motion\n"
        "\n"
        "This script will recursively search for valid python files\n"
        "If you only want to apply it to some file\n"
        "type: i-pi-style --file <file_path>\n"
        "example: i-pi-style -f ipi/engine/motion/motion.py\n",
    )

    parser.add_argument(
        "-p",
        "--path",
        type=str,
        default=None,
        help="Path to directory to look for python files recursively. Relative to the current path unless called with an absolute path.",
    )
    parser.add_argument(
        "-f",
        "--file_path",
        type=str,
        default=None,
        help="Filename on which enforce styling. Relative to the current path unless called with an absolute path. ",
    )
    parser.add_argument(
        "-c",
        "--check",
        action="store_true",
        help="Only checks style compliance but does not modify any file",
    )

    args = parser.parse_args()
    folder_path = args.path
    file_path = args.file_path
    check_mode = args.check

    main(check_mode, folder_path, file_path)
