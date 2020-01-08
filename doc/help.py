"""Help script which automatically generates help files.

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
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import time
from optparse import OptionParser

src_dir = "../"

sys.path.append(src_dir)

from ipi.inputs import barostats, cell, simulation, system, \
    ensembles, thermostats, forcefields, forces, atoms, beads, \
    prng, initializer, normalmodes, outputs, motion, smotion


time.sleep(1)
__all__ = ['help', 'objects']

objects = {'barostats': barostats.InputBaro(),
           'cell': cell.InputCell(),
           'simulation': simulation.InputSimulation(),
           'smotion': smotion.smotion.InputSmotion(),
           'remd': smotion.remd.InputReplicaExchange(),
           'meta': smotion.metad.InputMetaDyn(),
           'system': system.InputSystem(),
           'ensembles': ensembles.InputEnsemble(),
           'motion': motion.motion.InputMotion(),
           'dynamics': motion.dynamics.InputDynamics(),
           't_ramp': motion.ramp.InputTemperatureRamp(),
           'p_ramp': motion.ramp.InputPressureRamp(),
           'alchemy': motion.alchemy.InputAlchemy(),
           'planetary': motion.planetary.InputPlanetary(),
           'atomswap': motion.atomswap.InputAtomSwap(),
           'instanton': motion.instanton.InputInst(),
           'phonons': motion.phonons.InputDynMatrix(),
           'geop': motion.geop.InputGeop(),
           'neb': motion.neb.InputNEB(),
           'thermostats': thermostats.InputThermo(),
           'forcefield': forcefields.InputForceField(),
           'ffsocket': forcefields.InputFFSocket(),
           'fflj': forcefields.InputFFLennardJones(),
           'ffdebye': forcefields.InputFFDebye(),
           'ffquip': forcefields.InputFFQUIP(),
           'ffplumed': forcefields.InputFFPlumed(),
           'ffyaff': forcefields.InputFFYaff(),
           'ffsgdml': forcefields.InputFFsGDML(),
           'forcecomponent': forces.InputForceComponent(),
           'forces': forces.InputForces(),
           'atoms': atoms.InputAtoms(),
           'beads': beads.InputBeads(),
           'prng': prng.InputRandom(),
           'normalmodes': normalmodes.InputNormalModes(),
           'nmfrequencies': normalmodes.InputNMFrequencies(),
           'init_file': initializer.InputInitFile(),
           'init_pos': initializer.InputInitPositions(),
           'init_mom': initializer.InputInitMomenta(),
           'init_lab': initializer.InputInitLabels(),
           'init_mass': initializer.InputInitMasses(),
           'init_vel': initializer.InputInitVelocities(),
           'init_cell': initializer.InputInitCell(),
           'init_therm': initializer.InputInitThermo(),
           'initializer': initializer.InputInitializer(),
           'output': outputs.InputOutputs(),
           'properties': outputs.InputProperties(),
           'checkpoint': outputs.InputCheckpoint(),
           'trajectory': outputs.InputTrajectory()}

usage = "usage: python %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-x", action="store_true", dest="xml", default=False, help="write an xml help file")
parser.add_option("-l", action="store_true", dest="latex", default=False, help="write a latex help file")
parser.add_option("-n", action="store", type="int", dest="levels", help="number of levels depth to which data is printed out")
parser.add_option("-o", action="store", dest="prefix", help="Prefix for the output files", default="help")
parser.add_option("-i", action="store", dest="opt", help="Root object for the help files. Options: ['barostats', 'cell', 'simulation', 'ensembles', 'thermostats', 'interface', 'forces', 'atoms', 'beads', 'prng', 'output', 'trajectory', 'properties', 'checkpoint']", default='simulation')
parser.add_option("-r", action="store_true", dest="ref", default=False, help="add references to a latex help file. Ignored if -l is not present")
(options, args) = parser.parse_args()

if options.opt not in objects:
    raise ValueError("Option " + options.opt + " is not a viable tag name")


def help(latex=False, xml=False, levels=None, option='simulation', prefix="help", standalone=True):
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
        xml_output = open(prefix + ".xml", "w")
        xml_output.write(simrestart.help_xml(name=option, stop_level=levels))
    if latex:
        latex_output = open(prefix + ".tex", "w")
        latex_output.write(simrestart.help_latex(stop_level=levels, standalone=standalone))

if __name__ == '__main__':
    help(options.latex, options.xml, options.levels, options.opt, options.prefix, not options.ref)
