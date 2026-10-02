"""Checks the package version is exposed and matches its single source."""

import ipi
from ipi._version import __version__


def test_version_is_a_string():
    assert isinstance(__version__, str)
    assert __version__


def test_package_exposes_version():
    assert ipi.__version__ == __version__
