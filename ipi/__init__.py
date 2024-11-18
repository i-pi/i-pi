"""
The i-PI package.
"""

# expose some utility functions in a more direct way
import subprocess
from ipi.utils.parsing import read_output, read_trajectory
from ipi.utils.setup import install_driver

__all__ = [
    "clients",
    "engine",
    "inputs",
    "interfaces",
    "utils",
    "pes",
    "ipi_global_settings",
    "install_driver",
    "read_output",
    "read_trajectory",
]

ipi_global_settings = {"floatformat": "%16.8e"}

try:
    # Try to import the generated version file
    from ._version import __version__
except ImportError:
    # If the import fails, try to get the version from the git repository
    def get_git_version():
        try:
            # Run the 'git describe' command to get the current commit SHA
            result = subprocess.run(
                ["git", "describe", "--tags", "--always"],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                check=True,
            )
            return result.stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            return "unknown"  # Fallback if git is not available

    __version__ = get_git_version()

# Now you can use __version__ in your project
