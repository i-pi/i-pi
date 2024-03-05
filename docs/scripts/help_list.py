"""Help script which automatically generates help files.

This takes an output class specified by the user, and then uses the
automatic help generation functions to generate appropriate help files for this
class, giving information about the tags and the attributes of this class.

A full help message can be found by running 'python help.py -h' or
'python help.py --help'.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
from optparse import OptionParser
from ipi.engine.properties import Properties, Trajectories, help_latex, help_rst

src_dir = ".."

sys.path.append(src_dir)

__all__ = ["help_list", "list_objects"]

list_objects = {"property_list": Properties(), "trajectory_list": Trajectories()}

usage = "usage: python %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option(
    "-o",
    action="store",
    dest="prefix",
    help="Prefix for the output files",
    default="help",
)
parser.add_option(
    "-i",
    action="store",
    dest="opt",
    help="Root object for the help files. Options: ['property_list', 'trajectory_list']",
    default="property_list",
)
parser.add_option(
    "-r",
    action="store_true",
    dest="ref",
    default=False,
    help="If false, this creates a stand-alone document.",
)
(options, args) = parser.parse_args()

if options.opt not in list_objects:
    raise ValueError("Option " + options.opt + " is not a viable tag name")


def help_list(
    option="property_list", prefix="help", standalone=True, latex=True, rst=False
):
    """Writes the help file.

    Will write a latex file 'prefix.tex'. Can also print out
    sections of latex documents rather than entire
    documents, so that we can input them into other latex documents, such as
    the manual.

    Args:
       option: A string specifying which object will be used as the root object
          for the latex and xml files. Defaults to 'property_list'.
       prefix: File prefix for the output files. Defaults to 'help'.
       standalone: Boolean specifying whether the latex file will be a stand-alone
          document, or will instead be intended to be used in a larger document
          with cross references between the different objects.
    """

    simrestart = list_objects[option]
    if option == "property_list":
        idict = simrestart.property_dict
    elif option == "trajectory_list":
        idict = simrestart.traj_dict
    else:
        raise ValueError("Incorrect option specified.")

    if latex:
        latex_output = open(prefix + ".tex", "w")
        latex_output.write(help_latex(idict, standalone=standalone))
    if rst:
        rst_output = open(prefix + ".rst", "w")
        rst_output.write(help_rst(idict, standalone=standalone))


if __name__ == "__main__":
    help_list(options.opt, options.prefix, not options.ref)
