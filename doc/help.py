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


This takes an input class specified by the user, and then uses the
automatic help generation functions to generate appropriate help files for this
class, giving information about the tags and the attributes of this class.

There are several options that can be specified, including the depth of the tag
hierarchy that will be output, the output format and the output file name.

A full help message can be found by running 'python help.py -h' or
'python help.py --help'.

Note that any new input class type must be added to the objects 
dictionary and the latex help file must be added to the end of 
the manual.lyx file for it to be included in the automatic help generation.
If you do create a new input class type, please include this in the help string
for the -i option.

Functions:
   help: Writes the help file.
"""

import sys

src_dir = ".."

sys.path.append(src_dir)

from ipi.inputs import *
from ipi.utils.io.io_xml import *
from optparse import OptionParser

__all__ = ['help', 'objects']

objects = { 'barostats': barostats.InputBaro(), 
            'cell': cell.InputCell(),
            'simulation': simulation.InputSimulation(),
            'ensembles': ensembles.InputEnsemble(),
            'thermostats': thermostats.InputThermo(),
            'socket': forces.InputFBSocket(),
            'forces': forces.InputForces(),
            'atoms': atoms.InputAtoms(),
            'beads': beads.InputBeads(),
            'prng': prng.InputRandom(),
            'init_file': initializer.InputInitFile(),
            'init_pos': initializer.InputInitPositions(),
            'init_mom': initializer.InputInitMomenta(),
            'init_lab': initializer.InputInitLabels(),
            'init_mass': initializer.InputInitMasses(),
            'init_vel': initializer.InputInitVelocities(),
            'init_cell': initializer.InputInitCell(),
            'init_therm': initializer.InputInitThermo(),
            'initializer': initializer.InputInitializer(),
            'normal_modes': normalmodes.InputNormalModes(),
            'output': outputs.InputOutputs(),
            'properties': outputs.InputProperties(),
            'checkpoint': outputs.InputCheckpoint(),
            'trajectory': outputs.InputTrajectory() }

usage = "usage: python %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-x", action="store_true", dest = "xml", default=False, help="write an xml help file")
parser.add_option("-l", action="store_true", dest = "latex", default=False, help="write a latex help file")
parser.add_option("-n", action="store", type="int", dest="levels", help="number of levels depth to which data is printed out")
parser.add_option("-o", action="store", dest="prefix", help="Prefix for the output files", default="help")
parser.add_option("-i", action="store", dest="opt", help="Root object for the help files. Options: ['barostats', 'cell', 'simulation', 'ensembles', 'thermostats', 'interface', 'forces', 'atoms', 'beads', 'prng', 'output', 'trajectory', 'properties', 'checkpoint']", default='simulation')
parser.add_option("-r", action="store_true", dest = "ref", default=False, help="add references to a latex help file. Ignored if -l is not present")
(options, args) = parser.parse_args()

if options.opt not in objects:
   raise ValueError("Option " + options.opt + " is not a viable tag name")

def help(latex=False, xml=False, levels = None, option='simulation', prefix="help", standalone=True):
   """Writes the help file.

   Will write an xml file 'prefix.xml' if xml=True and a latex file 'prefix.tex'
   if latex=True. Will write out tags to a depth equal to the value of levels, 
   if it is specified, and will print using a root tag as specified by 
   the value of option. The output will be given by prefix.tex and/or 
   prefix.xml, if latex and/or xml is True respectively. Can also print out 
   sections of latex documents with cross-references rather than entire 
   documents, so that we can input them into other latex documents, such as
   the manual.

   Args:
      latex: Boolean specifying whether a latex file will be printed.
      xml: Boolean specifying whether an xml file will be printed.
      levels: An integer specifying how many layers of the hierarchy will be 
         printed. If not given, all layers will be printed.
      option: A string specifying which object will be used as the root object
         for the latex and xml files. Defaults to 'simulation'.
      prefix: File prefix for the output files. Defaults to 'help'.
      standalone: Boolean specifying whether the latex file will be a stand-alone
         document, or will instead be intended to be used in a larger document
         with cross references between the different objects.
   """

   simrestart = objects[option]
   
   if xml:
      xml_output = open(prefix + ".xml","w")
      xml_output.write(simrestart.help_xml(name=option, stop_level=levels))
   if latex:
      latex_output = open(prefix + ".tex","w")
      latex_output.write(simrestart.help_latex(stop_level=levels, standalone=standalone))

if __name__ == '__main__':
   help(options.latex, options.xml, options.levels, options.opt, options.prefix, not options.ref)
