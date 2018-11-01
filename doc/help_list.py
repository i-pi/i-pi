"""Help script which automatically generates help files.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.


This takes an output class specified by the user, and then uses the
automatic help generation functions to generate appropriate help files for this
class, giving information about the tags and the attributes of this class.

A full help message can be found by running 'python help.py -h' or
'python help.py --help'.

Functions:
   help: Writes the help file.
"""

import sys

src_dir = ".."

sys.path.append(src_dir)

from ipi.engine.properties import *
from ipi.utils.io.io_xml import *
from optparse import OptionParser

__all__ = ['help_list', 'list_objects']

list_objects = { 'property_list': Properties(), 
            'trajectory_list': Trajectories()}

usage = "usage: python %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-o", action="store", dest="prefix", help="Prefix for the output files", default="help")
parser.add_option("-i", action="store", dest="opt", help="Root object for the help files. Options: ['property_list', 'trajectory_list']", default='property_list')
parser.add_option("-r", action="store_true", dest = "ref", default=False, help="If false, this creates a stand-alone document.")
(options, args) = parser.parse_args()

if options.opt not in list_objects:
   raise ValueError("Option " + options.opt + " is not a viable tag name")

def help_list(option='property_list', prefix="help", standalone=True):
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
      idict = list_objects[option].property_dict
   elif option == "trajectory_list":
      idict = list_objects[option].traj_dict
   else:
      raise ValueError("Incorrect option specified.")
   
   latex_output = open(prefix + ".tex","w")
   latex_output.write(help_latex(idict, standalone=standalone))

if __name__ == '__main__':
   help_list(options.opt, options.prefix, not options.ref)
