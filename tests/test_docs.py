"""Ensure that docs can always be built."""

# This file is part of i-PI.
# i-PI Copyright (C) 2015 i-PI developers
# See the "licenses" directory for full license information.

import os
import subprocess
import shlex

import pytest


def run_command(cmd):
    """Runs `cmd` in doc directory."""
    cwd = os.getcwd()
    f_null = open(os.devnull, 'w')
    os.chdir(os.path.join(os.path.split(__file__)[0], "..", "doc"))
    ret = subprocess.call(shlex.split(cmd), stdout=f_null, stderr=f_null)
    f_null.close()
    os.chdir(cwd)
    return ret


@pytest.fixture
def distclean():
    """Prepare for documentation build testing."""
    run_command("make distclean")


@pytest.fixture
def clean():
    """Clean up the documentation build after testing."""
    yield
    run_command("make clean")


def test_make_aux(distclean):
    """doc: run make aux"""
    ret = run_command("make aux")
    assert ret == 0


def test_make(clean):
    """doc: run make"""
    ret = run_command("make")
    assert ret == 0
