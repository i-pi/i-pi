"""Classes to print info, warnings and errors to standard output during the simulation."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import traceback
import sys
import os

__all__ = [
    "Verbosity",
    "verbosity",
    "banner",
    "info",
    "warning",
    "get_identification_info",
]


VERB_QUIET = 0
VERB_LOW = 1
VERB_MEDIUM = 2
VERB_HIGH = 3
VERB_DEBUG = 4
VERB_TRACE = 5


class Verbosity(object):
    """Class used to determine what to print to standard output.

    Attributes:
        level: Determines what level of output to print.
    """

    lock = False
    level = VERB_LOW

    def __getattr__(self, name):
        """Determines whether a certain verbosity level is
        less than or greater than the stored value.

        Used to decide whether or not a certain info or warning string
        should be output.

        Args:
            name: The verbosity level at which the info/warning string
                will be output.
        """

        if name == "quiet":
            return self.level >= VERB_QUIET
        elif name == "low":
            return self.level >= VERB_LOW
        elif name == "medium":
            return self.level >= VERB_MEDIUM
        elif name == "high":
            return self.level >= VERB_HIGH
        elif name == "debug":
            return self.level >= VERB_DEBUG
        elif name == "trace":
            return self.level >= VERB_TRACE
        else:
            return super(Verbosity, self).__getattr__(name)

    def __setattr__(self, name, value):
        """Sets the verbosity level

        Args:
            name: The name of what to set. Should always be 'level'.
            value: The value to set the verbosity to.

        Raises:
            ValueError: Raised if either the name or the level is not
                a valid option.
        """

        if name == "level":
            if self.lock:
                # do not set the verbosity level if this is locked
                return
            if value == "quiet":
                level = VERB_QUIET
            elif value == "low":
                level = VERB_LOW
            elif value == "medium":
                level = VERB_MEDIUM
            elif value == "high":
                level = VERB_HIGH
            elif value == "debug":
                level = VERB_DEBUG
            elif value == "trace":
                level = VERB_TRACE
            else:
                raise ValueError(
                    "Invalid verbosity level " + str(value) + " specified."
                )
            super(Verbosity, self).__setattr__("level", level)
        else:
            super(Verbosity, self).__setattr__(name, value)


verbosity = Verbosity()


def read_git_file(filepath):
    """Reads and returns the content of a Git-related file."""
    try:
        with open(filepath, "r") as file:
            return file.read().strip()
    except FileNotFoundError:
        return None


def get_git_info():
    """
    Retrieves basic Git repository information by reading .git directory files.

    Returns
    -------
    dict or None
        A dictionary containing Git information or None if unable to read Git data.
    """
    base_path = os.path.abspath(os.path.join(__file__, "..", "..", "..")) + "/"
    git_dir = os.path.join(base_path, ".git")
    if not os.path.isdir(git_dir):
        return None

    # Read HEAD to get the current branch or commit hash
    head_content = read_git_file(os.path.join(git_dir, "HEAD"))
    if not head_content:
        return None

    branch_name = "unknown"
    last_commit = "unknown"
    remote_url = "unknown"
    commit_author = "unknown"
    commit_message = "unknown"

    # Parse the current branch name
    try:
        if head_content.startswith("ref:"):
            ref_path = os.path.join(git_dir, head_content.split(" ")[1])
            branch_name = os.path.basename(ref_path)
            last_commit = read_git_file(ref_path)
        else:
            # Detached HEAD state
            branch_name = "DETACHED"
            last_commit = head_content
    except:
        pass

    try:
        # Read the last commit details from logs if available
        logs_path = os.path.join(git_dir, "logs", "HEAD")
        # commit_author = commit_date = commit_message = None
        if os.path.exists(logs_path):
            with open(logs_path, "r") as logs_file:
                last_log = logs_file.readlines()[-1].split()
                commit_author = last_log[2]  # Simplified; may need adjustments
                commit_message = " ".join(last_log[7:])  # Adjust index as necessary
    except:
        pass

    try:
        # Get remote URL from the config file
        config_path = os.path.join(git_dir, "config")
        # remote_url = None
        if os.path.exists(config_path):
            with open(config_path, "r") as config_file:
                for line in config_file:
                    if line.strip().startswith("url ="):
                        remote_url = line.split("=")[1].strip()
                        break
    except:
        pass

    return {
        "branch_name": branch_name,
        "last_commit": last_commit,
        "remote_url": remote_url,
        "commit_author": commit_author,
        "commit_message": commit_message,
    }


def get_system_info():
    # Get the current working directory
    current_folder = os.getcwd()

    # Get the machine name (hostname)
    try:
        with open("/etc/hostname", "r") as file:
            machine_name = file.read().strip()
    except FileNotFoundError:
        machine_name = "Unknown"

    # Get the FQDN using a workaround if needed
    fqdn = machine_name  # Basic fallback, as direct FQDN might not be accessible

    # Get the operating system name and version from `/etc/os-release`
    os_name = "Unknown"
    os_version = "Unknown"
    try:
        with open("/etc/os-release", "r") as file:
            for line in file:
                if line.startswith("NAME="):
                    os_name = line.split("=")[1].strip().strip('"')
                elif line.startswith("VERSION="):
                    os_version = line.split("=")[1].strip().strip('"')
    except FileNotFoundError:
        pass

    # Get the processor name from `/proc/cpuinfo`
    processor = "Unknown"
    try:
        with open("/proc/cpuinfo", "r") as file:
            for line in file:
                if line.startswith("model name"):
                    processor = line.split(":")[1].strip()
                    break
    except FileNotFoundError:
        pass

    # Get the number of CPUs from `/proc/cpuinfo`
    num_nodes = 0
    try:
        with open("/proc/cpuinfo", "r") as file:
            num_nodes = sum(1 for line in file if line.startswith("processor"))
    except FileNotFoundError:
        pass

    # Get the user name from environment variables
    user_name = os.getenv("USER") or os.getenv("USERNAME") or "Unknown"

    return {
        "current_folder": current_folder,
        "machine_name": machine_name,
        "fqdn": fqdn,
        "os_name": os_name,
        "os_version": os_version,
        "processor": processor,
        "num_nodes": num_nodes,
        "user_name": user_name,
    }


def get_identification_info():
    """
    Collects and formats both Git and system information into a human-readable string.

    Returns
    -------
    str
        A formatted string with Git and system information.
    """
    git_info = get_git_info()
    system_info = get_system_info()

    info_string = ""

    if git_info:
        info_string += "# Git information:\n"
        info_string += f"#      Remote URL: {git_info['remote_url']:<24}\n"
        info_string += f"#          Branch: {git_info['branch_name']:<24}\n"
        info_string += f"#     Last Commit: {git_info['last_commit']:<24}\n"
        info_string += f"#   Commit Author: {git_info['commit_author'] or 'N/A':<24}\n"
        info_string += f"#  Commit Message: {git_info['commit_message'] or 'N/A':<24}\n"
    else:
        info_string += "# Unable to retrieve Git information.\n"

    info_string += "#\n"

    if system_info:
        info_string += "# System information:\n"
        info_string += f"#     Current Folder: {system_info['current_folder']}\n"
        info_string += f"#       Machine Name: {system_info['machine_name']}\n"
        info_string += f"#               FQDN: {system_info['fqdn']}\n"
        info_string += f"#   Operating System: {system_info['os_name']}\n"
        info_string += f"#         OS Version: {system_info['os_version']}\n"
        info_string += f"#          Processor: {system_info['processor']}\n"
        info_string += f"#     Number of CPUs: {system_info['num_nodes']}\n"
        info_string += f"#          User Name: {system_info['user_name']}\n"
    else:
        info_string += "# Unable to retrieve system information.\n"

    return info_string


def banner():
    """Prints out a banner."""

    print(
        r"""
 ____       ____       ____       ____
/    \     /    \     /    \     /    \
|  #################################  |
\__#_/     \____/     \____/     \_#__/
   #    _        _______  _____    #
   #   (_)      |_   __ \|_   _|   #      -*-       v 3.0      -*-
   #   __  ______ | |__) | | |     #
   Y  [  ||______||  ___/  | |     #      A Universal Force Engine
  0 0  | |       _| |_    _| |_    #
   #  [___]     |_____|  |_____|   #
 __#_       ____       ____       _#__
/  # \     /    \     /    \     / #  \
|  #################################  |
\____/     \____/     \____/     \____/

    """
    )

    info_string = get_identification_info()
    print(info_string)


def info(text="", show=True):
    """Prints a message.

    Args:
        text: The text of the information message.
        show: A boolean describing whether or not the message should be
            printed.
    """

    if not show:
        return
    print(text)


def warning(text="", show=True):
    """Prints a warning message.

    Same as info, but with a "!W!" prefix and optionally printing a stack trace.

    Args:
        text: The text of the information message.
        show: A boolean describing whether or not the message should be
            printed.
    """

    if not show:
        return
    if verbosity.trace:
        traceback.print_stack(file=sys.stdout)
    print((" !W! " + text))
