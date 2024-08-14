"""Classes to print info, warnings and errors to standard output during the simulation."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import traceback
import sys

__all__ = ["Verbosity", "verbosity", "banner", "info", "warning"]


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
