"""Tests for ipi version handling."""

import pytest
import sys
from importlib.metadata import PackageNotFoundError


def test_version_from_package():
    """Test that version is retrieved from installed package."""
    import ipi

    # When installed via pip, version should be accessible
    assert ipi.__version__ is not None
    assert len(ipi.__version__) > 0
    assert ipi.__version__ != "unknown"


def test_version_fallback_when_not_installed(monkeypatch):
    """Test fallback to hardcoded version when PackageNotFoundError is raised.
    
    This test verifies that when importlib.metadata.version() fails,
    the module still has a valid version (not "unknown").
    """

    def mock_version(package_name):
        """Mock version() to simulate package not being installed."""
        if package_name == "ipi":
            raise PackageNotFoundError("ipi")
        return package_name

    monkeypatch.setattr("importlib.metadata.version", mock_version)

    # Remove ipi from sys.modules to force reimport with the mocked version()
    if "ipi" in sys.modules:
        del sys.modules["ipi"]

    import ipi as ipi_reimported

    # The key assertion: version should NOT be "unknown"
    # This proves the fallback to hardcoded version is working
    assert ipi_reimported.__version__ != "unknown"
    assert len(ipi_reimported.__version__) > 0
