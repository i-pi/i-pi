#!/usr/bin/env python3

"""
This script cleans python files in-place by calling autopep8 tool.
With option --path, it searches for python files
recursively in the given directory. With option --files
it applies autopep8 to listed python files.
Prior to running the script, autopep8 must be installed.

autopep8 is called with the following options:
-i (in-place, modifies files in-place)
-v (verbose)
--select errors (only selected pep8 issues are fixed)
--pep8-passes 2000 (avoids falling into infinite loop by autopep8).

In the check mode, the tool uses pycodestyle API.
http://pycodestyle.readthedocs.io/en/latest/#

For a full list of issues which can be fixed by autopep8, consult:
https://pypi.python.org/pypi/autopep8
"""

import argparse
import os
import subprocess
import re
import sys
import autopep8
import pycodestyle


if __name__ == "__main__":
    # Description of the program after using -h
    parser = argparse.ArgumentParser(
        description="Pepper executes autopep8 "
        "cleaning of python files according to pep8 style guide. "
        "To run it autopep8 must be installed: "
        "https://pypi.python.org/pypi/autopep8  "
        "With --path option, it recursively looks "
        "for all python files in a given directory, "
        "so for example, if you want to clean your i-pi repository, "
        "type: pepper --path your_ipi_root . "
        "Pepper with recursively search for valid python files and "
        "will clean them. "
        "If you only want to apply it to some files, use --files "
        "option, pepper will check if they are python files!"
    )
    # There is only one positional argument
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-p",
        "--path",
        type=str,
        metavar="PATH",
        help="Path to directory "
        "where pepper will look for python files recursively. "
        "May not be used with --files",
    )
    group.add_argument(
        "-f",
        "--files",
        type=str,
        nargs="+",
        help="list od files on which pepper will execute "
        "autopep8. May not be used with --path",
    )
    parser.add_argument(
        "-v",
        "--verbosity",
        choices=["silent", "low", "medium", "high"],
        default="medium",
        help="sets level of verbosity. silent will not print "
        "anything. low prints messages from this script "
        "and no autopep8 output, "
        "medium prints filenames checked by pepper"
        "and informs of any unfixed problems, "
        "high prints everything from autopep8 output",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="checks pep8 compliance, does not modify files",
    )
    args = parser.parse_args()
    path = args.path
    files = args.files
    verbosity = args.verbosity
    is_in_check_mode = args.check

    # General arguments for pep8
    styles_to_be_corrected = [
        "E101",
        "E11",
        "E121",
        "E122",
        "E123",
        "E124",
        "E125",
        "E126",
        "E127",
        "E128",
        "E20",
        "E211",
        "E22",
        "E224",
        "E226",
        "E227",
        "E228",
        "E231",
        "E241",
        "E242",
        "E251",
        "E26",
        "E265",
        "E27",
        "E301",
        "E302",
        "E303",
        "E304",
        "E306",
        "E401",
        "E502",
        "W291",
        "W292",
        "W29",
        "W391",
    ]
    # Must be written as string since it is passed to command
    number_of_pep8_passes = "2000"
    # in-place, verbose, select only styles which are listed above,
    # pep8-passes limit is to avoid infinite loop
    autopep8_args = [
        "-i",
        "-v",
        "--select",
        ",".join(styles_to_be_corrected),
        "--pep8-passes",
        number_of_pep8_passes,
    ]

    if files is not None:
        python_files = []
        for filename in files:
            if autopep8.is_python_file(filename):
                python_files = python_files + [filename]
            else:
                if verbosity != "silent":
                    print(filename, "is not a python file, skipping ")
        # execute autopep8 only on python files, without being recursive
        autopep8_args = autopep8_args + python_files
        if not python_files:
            print("No python files to process.")
            sys.exit()
        if verbosity != "silent":
            print("Running autopep8 on the files: ", " ".join(python_files))
    else:
        # perform recursive search in the given directory
        if os.path.isdir(path):
            os.chdir(path)
            autopep8_args = autopep8_args + ["-r"] + ["."]
            if verbosity != "silent":
                print("Running autopep8 recursively in the directory: ", os.getcwd())
        else:
            print("The given directory does not exist: ", path)
            sys.exit()

    if is_in_check_mode:
        # perform check, do not change files
        print("Performing check of PEP-8 compliance")
        style_check = pycodestyle.StyleGuide(select=styles_to_be_corrected)
        if files is not None:
            for filename in files:
                style_check.input_file(filename)
        else:
            style_check.input_dir(".")
    else:
        # use autopep8 to change files in-place
        process = subprocess.Popen(
            ["autopep8"] + autopep8_args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        # Line buffer to save the last unfixed issue: autopep8 prints unfixed issues at
        # consecutive iterations on a given file. This generates a lot of messages.
        # We are only interested in the last message and print it if number of issues exceeds 0.
        # This means that some unfixed issues are left in the processed file.
        issue_line_buffer = ""
        while process.poll() is None:
            # We must strip, otherwise we get double newline
            line = process.stdout.readline().rstrip()
            if re.match(r"\[Errno.*\]", line):
                print(line)
            else:
                if verbosity == "high":
                    print(line)
                elif verbosity == "medium":
                    if re.match(r"\[file:.*\]", line):
                        # new file is processed. If something left unfixed, print message now
                        if issue_line_buffer != "":
                            print(issue_line_buffer)
                        # Pattern: [file:filename]
                        filename_line = re.search(r"\[file:(\S+)\]", line)
                        # Print only filename
                        print(filename_line.group(1))
                    elif re.match(".*0 issue.*", line):
                        # there are 0 issues left, clean the buffer
                        issue_line_buffer = ""
                    elif re.match(".*[1-9] issue.*", line):
                        # number of issues is greater than 0
                        issue_line = line
                        # Reformat this line:
                        # --->  6 issue(s) to fix {'E266': set([131, 931]), 'E265': set([512, 164, 518, 519])}
                        # to achieve that message and formatting:
                        # --->  6 issue(s) not fixed in file ./regtest.py
                        # 'E266' at lines: 131, 931
                        # 'E265' at lines: 512, 164, 518, 519
                        issue_line = issue_line.replace(
                            "to fix", "not fixed in file " + filename_line.group(1)
                        )
                        issue_line = issue_line.replace("{", "\n      ")
                        issue_line = issue_line.replace("]),", "\n     ")
                        issue_line = issue_line.replace(": set([", " at lines: ")
                        issue_line = issue_line.replace("])", "")
                        issue_line = issue_line.replace("}", "")
                        issue_line_buffer = issue_line
                # if verbosity is silent or low, do not print output from autopep8
        if verbosity == "medium":
            # end of processing. If something left unfixed, print the message now
            if issue_line_buffer == "":
                print(issue_line_buffer)
    if verbosity != "silent":
        print("autopep8 terminated")
