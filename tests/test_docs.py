"""Ensure that docs can always be built."""

# This file is part of i-PI.
# i-PI Copyright (C) 2015 i-PI developers
# See the "licenses" directory for full license information.

import os
import subprocess
import shlex

from nose import with_setup

FNULL = open(os.devnull, 'w')


def run_command(cmd):
    """Runs `cmd` in doc directory."""
    cwd = os.getcwd()
    os.chdir(os.path.join(os.path.split(__file__)[0], "..", "doc"))
    ret = subprocess.call(shlex.split(cmd), stdout=FNULL, stderr=FNULL)
    os.chdir(cwd)
    return ret


def distclean():
    """Prepare for documentation build testing."""
    run_command("make distclean")


def clean():
    """Clean up the documentation build after testing."""
    run_command("make clean")


@with_setup(distclean, None)
def test_make_aux():
    """doc: run make aux"""
    ret = run_command("make aux")
    assert ret == 0


@with_setup(None, clean)
def test_make():
    """doc: run make"""
    ret = run_command("make")
    assert ret == 0
