"""Holds the class which computes important properties of the system, and
prepares them for output.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.messages import verbosity, info, warning
from ipi.utils.depend import *
from ipi.utils.units import Constants, unit_to_internal
from ipi.utils.mathtools import logsumlog, h2abc_deg
from ipi.utils.io.inputs import io_xml
from ipi.engine.atoms import *
from ipi.engine.cell import *
from ipi.engine.ensembles import *
from ipi.engine.forces import *


__all__ = ["Properties", "Trajectories", "getkey", "getall", "help_latex"]


def getkey(pstring):
    """Strips units and argument lists from a property/trajectory keyword.

    Args:
       pstring: The string input by the user that specifies an output,
          which in general will specify units and argument lists.

    Returns: A string giving the keyword for the property, stripped of the
       argument lists and units key words.
    """

    pa = pstring.find("(")
    if pa < 0:
        pa = len(pstring)
    pu = pstring.find("{")
    if pu < 0:
        pu = len(pstring)
    return pstring[0 : min(pa, pu)].strip()


def getall(pstring):
    """Returns the keyword, units and argument list separately.

    Args:
       pstring: The string input by the user that specifies an output,
          which in general will specify units and argument lists.

    Returns: A tuple giving the keyword for the property, and its units
       argument list and key word argument list.
    """

    unit = ""
    arglist = ()
    kwarglist = {}
    unstart = len(pstring)
    argstart = unstart

    if "}" in pstring:
        # the property has a user-defined unit
        unstart = pstring.find("{")
        unstop = pstring.find("}", unstart)
        if unstop == -1:
            raise ValueError("Incorrect format in units specification " + pstring)
        unit = pstring[unstart + 1 : unstop]
    if "(" in pstring:
        # If the property has additional arguments
        argstart = pstring.find("(")
        argstop = pstring.find(")", argstart)
        if argstop == -1:
            raise ValueError("Incorrect format in argument list " + pstring)

        argstr = pstring[argstart : argstop + 1]
        arglist = io_xml.read_tuple(argstr, delims="()", split=";", arg_type=str)
        for arg in arglist:
            # If a keyword argument is used
            equals = arg.find("=")
            if equals >= 0:
                kwarglist[arg[0:equals].strip()] = arg[equals + 1 :].strip()
                arglist = tuple(a for a in arglist if not a == arg)

    pstring = pstring[
        0 : min(unstart, argstart)
    ].strip()  # strips the arguments from pstring name

    return (pstring, unit, arglist, kwarglist)


def help_latex(idict, standalone=True):
    """Function to generate a LaTeX formatted string.

    Can be used in the manual to list the different available outputs.

    Args:
       idict: Either property_dict or traj_dict, to be used to
          generate the help file.
       standalone: A boolean giving whether the latex file produced will be a
          stand-alone document, or will be intended as a section of a larger
          document with cross-references between the different sections.

    Returns:
       A LaTeX formatted string.
    """

    rstr = ""
    if standalone:
        # assumes that it is a stand-alone document, so must have document
        # options.
        rstr += r"\documentclass[12pt,fleqn]{report}"
        rstr += r"""
\usepackage{etoolbox}
\usepackage{suffix}

\newcommand{\ipiitem}[3]{%
\ifblank{#1}{}{\ifstrequal{#1}{\underline{}}{}{
{\noindent\textbf{#1}:\rule{0.0pt}{1.05\baselineskip}\quad}}}% uses a strut to add a bit of vertical space
{#2}\parskip=0pt\par
\ifblank{#3}{}%
{ {\hfill\raggedleft\textit{\small #3}\par} }
}
"""
        rstr += "\n" + r"\begin{document}\n"
        rstr += "The following are the different allowable ouputs:\n" + r"\par"

    for out in sorted(idict):
        rstr += r"\ipiitem{" + out + "}"
        if "longhelp" in idict[out]:
            rstr += "{" + idict[out]["longhelp"] + "}"
        else:
            rstr += "{" + idict[out]["help"] + "}"

        # see if there are additional attributes to print out
        xstr = ""
        if (
            "dimension" in idict[out] and idict[out]["dimension"] != "undefined"
        ):  # doesn't print out dimension if not necessary.
            xstr += "dimension: " + idict[out]["dimension"] + "; "
        if "size" in idict[out]:
            xstr += "size: " + str(idict[out]["size"]) + "; "
        rstr += "{" + xstr + "}"

    if standalone:
        # ends the created document if it is not part of a larger document
        rstr += r"\end{document}"

    # Some escape characters are necessary for the proper latex formatting
    rstr = rstr.replace("_", "\\_")
    rstr = rstr.replace("\\\\_", "\\_")
    rstr = rstr.replace("...", "\\ldots ")
    rstr = rstr.replace("<", "$<$")
    rstr = rstr.replace(">", "$>$")
    rstr = rstr.replace("[", "$[$")
    rstr = rstr.replace("]", "$]$")

    return rstr


class Properties(dobject):

    """A proxy to compute and output properties of the system.

    Takes the fundamental properties calculated during the simulation, and
    prepares them for output. It also contains simple algorithms to calculate
    other properties not calculated during the simulation itself, so that
    these can also be output.

    Attributes:
       fd_delta: A float giving the size of the finite difference
          parameter used in the Yamamoto kinetic energy estimator. Defaults
          to _DEFAULT_FINDIFF.
       _DEFAULT_FDERROR: A float giving the size of the minimum precision
          allowed for the finite difference calculation in the Yamamoto kinetic
          energy estimator.
       _DEFAULT_MINFID: A float giving the maximum displacement in the Yamamoto
          kinetic energy estimator.
       dbeads: A dummy Beads object used in the Yamamoto kinetic energy
          estimator.
       dforces: A dummy Forces object used in the Yamamoto kinetic energy
          estimator.
       system: The System object containing the data to be output.
       ensemble: An ensemble object giving the objects necessary for producing
          the correct ensemble.
       beads: A beads object giving the atoms positions.
       nm: A normal modes object giving the normal mode representation.
       cell: A cell object giving the system box.
       forces: A forcefield object giving the force calculator for each
          replica of the system.
       property_dict: A dictionary containing all the properties that can be
          output.
    """

    _DEFAULT_FINDIFF = 1e-4
    _DEFAULT_FDERROR = 1e-6
    _DEFAULT_MINFID = 1e-7

    def __init__(self):
        """Initialises Properties."""

        self.property_dict = {
            "step": {
                "dimension": "number",
                "help": "The current simulation time step.",
                "func": (lambda: (1 + self.simul.step)),
            },
            "time": {
                "dimension": "time",
                "help": "The elapsed simulation time.",
                "func": (lambda: self.ensemble.time),
            },
            "temperature": {
                "dimension": "temperature",
                "help": "The current temperature, as obtained from the MD kinetic energy.",
                "longhelp": """The current temperature, as obtained from the MD kinetic energy of the (extended)
                                      ring polymer. Takes optional arguments 'atom', 'bead' or 'nm'.  'atom' can be either an
                                      atom label or an index (zero-based) to specify which species or individual atom
                                      to output the temperature of. If not specified, all atoms are used and averaged.
                                      'bead' or 'nm' specify whether the temperature should be computed for a single bead
                                      or normal mode.""",
                "func": self.get_temp,
            },
            "density": {
                "dimension": "density",
                "help": "The mass density of the physical system.",
                "func": (lambda: self.beads.m.sum() / self.cell.V),
            },
            "volume": {
                "dimension": "volume",
                "help": "The volume of the cell box.",
                "func": (lambda: self.cell.V),
            },
            "cell_h": {
                "dimension": "length",
                "help": "The simulation cell as a matrix. Returns the 6 non-zero components in the form [xx, yy, zz, xy, xz, yz].",
                "size": 6,
                "func": (lambda: self.tensor2vec(self.cell.h)),
            },
            "cell_abcABC": {
                "dimension": "undefined",
                "help": "The lengths of the cell vectors and the angles between them in degrees as a list of the form [a, b, c, A, B, C]",
                "longhelp": """The lengths of the cell vectors and the angles between them in degrees as a list of the
                      form [a, b, c, A, B, C], where A is the angle between the sides of length b and c in degrees, and B and C
                      are defined similarly. Since the output mixes different units, a, b and c can only be output in bohr.""",
                "size": 6,
                "func": (lambda: np.asarray(h2abc_deg(self.cell.h))),
            },
            "conserved": {
                "dimension": "energy",
                "help": "The value of the conserved energy quantity per bead.",
                "func": (lambda: self.ensemble.econs / float(self.beads.nbeads)),
            },
            "ensemble_lp": {
                "dimension": "undefined",
                "help": "The log of the ensemble probability",
                "func": (lambda: self.ensemble.lpens),
            },
            "ensemble_temperature": {
                "dimension": "temperature",
                "help": "The target temperature for the current ensemble",
                "func": (lambda: self.ensemble.temp),
            },
            "ensemble_pressure": {
                "dimension": "pressure",
                "help": "The target pressure for the current ensemble",
                "func": (lambda: self.ensemble.pext),
            },
            "hweights_component": {
                "dimension": "",
                "help": "The weight associated to the one part of the hamiltonian. ",
                "longhelp": """The weight associated one part of the hamiltonian. Takes one mandatory
                         argument index (zero-based) that indicates for which component of the hamiltonian the weight must be returned. """,
                "func": (lambda index: self.ensemble.hweights[int(index)]),
            },
            "ensemble_bias": {
                "dimension": "energy",
                "help": "The bias applied to the current ensemble",
                "func": (lambda: self.ensemble.bias.pot / self.beads.nbeads),
            },
            "bweights_component": {
                "dimension": "",
                "help": "The weight associated to the one part of the hamiltonian. ",
                "longhelp": """The weight associated one part of the hamiltonian. Takes one mandatory
                         argument index (zero-based) that indicates for which component of the hamiltonian the weight must be returned. """,
                "func": (lambda index: self.ensemble.bweights[int(index)]),
            },
            #      "ensemble_logweight":  {  "dimension": "",
            #                       "help" : "The (log) weight of the configuration in the biassed ensemble",
            #                       "func": (lambda: self.ensemble.bias/(Constants.kb*self.ensemble.temp)) },
            "potential": {
                "dimension": "energy",
                "help": "The physical system potential energy.",
                "longhelp": """The physical system potential energy. With the optional argument 'bead'
                         will print the potential associated with the specified bead.""",
                "func": (
                    lambda bead="-1": self.forces.pot / self.beads.nbeads
                    if int(bead) < 0
                    else self.forces.pots[int(bead)]
                ),
            },
            "potential_opsc": {
                "dimension": "energy",
                "help": "The Suzuki-Chin operator estimator for the potential energy of the physical system.",
                "func": (
                    lambda: 2.0 / self.beads.nbeads * np.sum(self.forces.pots[::2])
                ),
            },
            "potential_tdsc": {
                "dimension": "energy",
                "help": "The Suzuki-chin thermodyanmic estimator for the potential energy of the physical system.",
                "func": self.get_scpottd,
            },
            "pot_component": {
                "dimension": "energy",
                "help": "The contribution to the system potential from one of the force components. ",
                "longhelp": """The contribution to the system potential from one of the force components.
                       Takes one mandatory argument index (zero-based) that indicates which component of the
                       potential must be returned. The optional argument 'bead' will print the potential associated
                       with the specified bead. If the potential is weighed, the weight will be applied. """,
                "func": (
                    lambda index, bead="-1": self.forces.pots_component(
                        int(index)
                    ).sum()
                    / self.beads.nbeads
                    if int(bead) < 0
                    else self.forces.pots_component(int(index))[int(bead)]
                ),
            },
            "pot_component_raw": {
                "dimension": "energy",
                "help": "The contribution to the system potential from one of the force components. ",
                "longhelp": """The contribution to the system potential from one of the
                       force components. Takes one mandatory argument index (zero-based) that indicates
                       which component of the potential must be returned. The optional argument 'bead'
                       will print the potential associated with the specified bead. Potential weights
                       will not be applied. """,
                "func": (
                    lambda index, bead="-1": self.forces.pots_component(
                        int(index), False
                    ).sum()
                    / self.beads.nbeads
                    if int(bead) < 0
                    else self.forces.pots_component(int(index), False)[int(bead)]
                ),
            },
            "forcemod": {
                "dimension": "force",
                "help": "The modulus of the force.",
                "longhelp": """The modulus of the force. With the optional argument 'bead'
                       will print the force associated with the specified bead.""",
                "func": (
                    lambda bead="-1": np.linalg.norm(self.forces.f) / self.beads.nbeads
                    if int(bead) < 0
                    else np.linalg.norm(self.forces.f[int(bead)])
                ),
            },
            "spring": {
                "dimension": "energy",
                "help": "The total spring potential energy between the beads of all the ring polymers in the system.",
                "func": (lambda: self.nm.vspring / self.beads.nbeads),
            },
            "kinetic_md": {
                "dimension": "energy",
                "help": "The kinetic energy of the (extended) classical system.",
                "longhelp": """The kinetic energy of the (extended) classical system.
                       Takes optional arguments 'atom', 'bead' or 'nm'.  'atom' can be either an
                       atom label or an index (zero-based) to specify which species or individual atom
                       to output the kinetic energy of. If not specified, all atoms are used and averaged.
                       'bead' or 'nm' specify whether the kinetic energy should be computed for a single bead
                       or normal mode. If not specified, all atoms/beads/nm are used.""",
                "func": self.get_kinmd,
            },
            "kinetic_cv": {
                "dimension": "energy",
                "help": "The centroid-virial quantum kinetic energy of the physical system.",
                "longhelp": """The centroid-virial quantum kinetic energy of the physical system.
                      Takes an argument 'atom', which can be either an atom label or index (zero based)
                      to specify which species to find the kinetic energy of. If not specified, all atoms are used.""",
                "func": self.get_kincv,
            },
            "kinetic_td": {
                "dimension": "energy",
                "help": "The primitive quantum kinetic energy of the physical system.",
                "longhelp": """The primitive quantum kinetic energy of the physical system.
                      Takes an argument 'atom', which can be either an atom label or index (zero based)
                      to specify which species to find the kinetic energy of. If not specified, all atoms are used.""",
                "func": self.get_kintd,
            },
            "kinetic_prsc": {
                "dimension": "energy",
                "help": "The Suzuki-Chin primitive estimator of the quantum kinetic energy of the physical system",
                "func": self.get_sckinpr,
            },
            "kinetic_tdsc": {
                "dimension": "energy",
                "help": "The Suzuki-Chin centroid-virial thermodynamic estimator of the quantum kinetic energy of the physical system.",
                "longhelp": """The Suzuki-Chin centroid-virial thermodynamic estimator of the quantum
                      kinetic energy of the physical system. Takes an argument 'atom', which can be either
                      an atom label or index (zero based) to specify which species to find the kinetic energy
                      of. If not specified, all atoms are used.""",
                "func": self.get_sckintd,
            },
            "kinetic_opsc": {
                "dimension": "energy",
                "help": "The Suzuki-Chin centroid-virial operator estimator of the quantum kinetic energy of the physical system.",
                "longhelp": """The centroid-virial quantum kinetic energy of the physical system.
                      Takes an argument 'atom', which can be either an atom label or index (zero based)
                      to specify which species to find the kinetic energy of. If not specified, all atoms are used.""",
                "func": self.get_sckinop,
            },
            "kinetic_tens": {
                "dimension": "energy",
                "help": "The centroid-virial quantum kinetic energy tensor of the physical system.",
                "longhelp": """The centroid-virial quantum kinetic energy tensor of the physical system.
                      Returns the 6 independent components in the form [xx, yy, zz, xy, xz, yz]. Takes an
                      argument 'atom', which can be either an atom label or index (zero based) to specify
                      which species to find the kinetic tensor components of. If not specified, all atoms are used.""",
                "size": 6,
                "func": self.get_ktens,
            },
            "kinetic_ij": {
                "dimension": "energy",
                "help": "The centroid-virial off-diagonal quantum kinetic energy tensor of the physical system.",
                "longhelp": """The centroid-virial off-diagonal quantum kinetic energy tensor of the physical system.
                      This computes the cross terms between atoms i and atom j, whose average is  <p_i*p_j/(2*sqrt(m_i*m_j))>.
                      Returns the 6 independent components in the form [xx, yy, zz, xy, xz, yz]. Takes arguments 'i' and 'j',
                       which give the indices of the two desired atoms.""",
                "size": 6,
                "func": self.get_kij,
            },
            "r_gyration": {
                "dimension": "length",
                "help": "The average radius of gyration of the selected ring polymers.",
                "longhelp": """The average radius of gyration of the selected ring polymers. Takes an
                      argument 'atom', which can be either an atom label or index (zero based) to specify which
                      species to find the radius of gyration of. If not specified, all atoms are used and averaged.""",
                "func": self.get_rg,
            },
            "atom_x": {
                "dimension": "length",
                "help": "The position (x,y,z) of a particle given its index.",
                "longhelp": """The position (x,y,z) of a particle given its index. Takes arguments index
                       and bead (both zero based). If bead is not specified, refers to the centroid.""",
                "size": 3,
                "func": (
                    lambda atom="", bead="-1": self.get_atom_vec(
                        self.beads.q, atom=atom, bead=bead
                    )
                ),
            },
            "atom_x_path": {
                "dimension": "length",
                "help": "The positions of all the beads of a particle given its index.",
                "longhelp": """The positions of all the beads of a particle given its index. Takes arguments index
                       and bead (both zero based). If bead is not specified, refers to the centroid.""",
                "func": (
                    lambda atom="": np.asarray(
                        [
                            self.get_atom_vec(self.beads.q, atom=atom, bead=j)
                            for j in range(self.beads.nbeads)
                        ]
                    ).flatten()
                ),
            },
            "atom_f_path": {
                "dimension": "length",
                "help": "The forces acting on all the beads of a particle given its index.",
                "longhelp": """The forces acting on all the beads of a particle given its index. Takes arguments index
                       and bead (both zero based). If bead is not specified, refers to the centroid.""",
                "func": (
                    lambda atom="": np.asarray(
                        [
                            self.get_atom_vec(self.forces.f, atom=atom, bead=j)
                            for j in range(self.beads.nbeads)
                        ]
                    ).flatten()
                ),
            },
            "atom_v": {
                "dimension": "velocity",
                "help": "The velocity (x,y,z) of a particle given its index.",
                "longhelp": """The velocity (x,y,z) of a particle given its index. Takes arguments index
                       and bead (both zero based). If bead is not specified, refers to the centroid.""",
                "size": 3,
                "func": (
                    lambda atom="", bead="-1": self.get_atom_vec(
                        self.beads.p / self.beads.m3, atom=atom, bead=bead
                    )
                ),
            },
            "vcom": {
                "dimension": "velocity",
                "help": "The COM velocity (x,y,z) of the system or a chosen species.",
                "longhelp": """The center of mass velocity (x,y,z) of the system or of a species. Takes arguments label
                       (default to all species) and bead (zero based). If bead is not specified, refers to the centroid.""",
                "size": 3,
                "func": self.get_vcom,
            },
            "atom_p": {
                "dimension": "momentum",
                "help": "The momentum (x,y,z) of a particle given its index.",
                "longhelp": """The momentum (x,y,z) of a particle given its index. Takes arguments index
                      and bead (both zero based). If bead is not specified, refers to the centroid.""",
                "size": 3,
                "func": (
                    lambda atom="", bead="-1": self.get_atom_vec(
                        self.beads.p, atom=atom, bead=bead
                    )
                ),
            },
            "atom_f": {
                "dimension": "force",
                "help": "The force (x,y,z) acting on a particle given its index.",
                "longhelp": """The force (x,y,z) acting on a particle given its index. Takes arguments index
                      and bead (both zero based). If bead is not specified, refers to the centroid.""",
                "size": 3,
                "func": (
                    lambda atom="", bead="-1": self.get_atom_vec(
                        self.forces.f, atom=atom, bead=bead
                    )
                ),
            },
            "stress_md": {
                "dimension": "pressure",
                "size": 6,
                "help": "The total stress tensor of the (extended) classical system.",
                "longhelp": """The total stress tensor of the (extended) classical system. Returns the 6
                      independent components in the form [xx, yy, zz, xy, xz, yz].""",
                "func": (
                    lambda: self.tensor2vec(
                        (self.forces.vir + self.nm.kstress) / self.cell.V
                    )
                ),
            },
            "pressure_md": {
                "dimension": "pressure",
                "help": "The pressure of the (extended) classical system.",
                "func": (
                    lambda: np.trace(
                        (self.forces.vir + self.nm.kstress) / (3.0 * self.cell.V)
                    )
                ),
            },
            "kstress_md": {
                "dimension": "pressure",
                "size": 6,
                "help": "The kinetic stress tensor of the (extended) classical system.",
                "longhelp": """The kinetic stress tensor of the (extended) classical system. Returns the 6
                      independent components in the form [xx, yy, zz, xy, xz, yz].""",
                "func": (lambda: self.tensor2vec(self.nm.kstress / self.cell.V)),
            },
            "virial_fq": {
                "dimension": "energy",
                "size": 1,
                "help": "The scalar product of force and position.",
                "longhelp": """Returns the scalar product of force and positions. Useful to compensate for
                          the harmonic component of a potential. Gets one argument 'ref' that should be a filename for a
                          reference configuration, in the style of the FFDebye geometry input, and one that contains the input units.""",
                "func": self.get_fqvirial,
            },
            "virial_md": {
                "dimension": "pressure",
                "size": 6,
                "help": "The virial tensor of the (extended) classical system.",
                "longhelp": """The virial tensor of the (extended) classical system. Returns the 6
                      independent components in the form [xx, yy, zz, xy, xz, yz].""",
                "func": (lambda: self.tensor2vec(self.forces.vir / self.cell.V)),
            },
            "stress_cv": {
                "dimension": "pressure",
                "size": 6,
                "help": "The total quantum estimator for the stress tensor of the physical system.",
                "longhelp": """The total quantum estimator for the stress tensor of the physical system. Returns the
                      6 independent components in the form [xx, yy, zz, xy, xz, yz].""",
                "func": (
                    lambda: self.tensor2vec(self.forces.vir + self.kstress_cv())
                    / (self.cell.V * self.beads.nbeads)
                ),
            },
            "pressure_tdsc": {
                "dimension": "pressure",
                "help": "The Suzuki-Chin thermodynamic estimator for pressure of the physical system.",
                "func": (
                    lambda: np.trace(
                        self.forces.vir + self.forces.virsc + self.kstress_sctd()
                    )
                    / (3.0 * self.cell.V * self.beads.nbeads)
                ),
            },
            "vir_tdsc": {
                "dimension": "pressure",
                "help": "The Suzuki-Chin thermodynamic estimator for pressure of the physical system.",
                "func": (
                    lambda: np.trace(self.forces.vir + self.forces.virsc)
                    / (3.0 * self.cell.V * self.beads.nbeads)
                ),
            },
            "kstress_tdsc": {
                "dimension": "pressure",
                "help": "The Suzuki-Chin thermodynamic estimator for pressure of the physical system.",
                "func": (
                    lambda: np.trace(self.kstress_sctd())
                    / (3.0 * self.cell.V * self.beads.nbeads)
                ),
            },
            "pressure_cv": {
                "dimension": "pressure",
                "help": "The quantum estimator for pressure of the physical system.",
                "func": (
                    lambda: np.trace(self.forces.vir + self.kstress_cv())
                    / (3.0 * self.cell.V * self.beads.nbeads)
                ),
            },
            "kstress_cv": {
                "dimension": "pressure",
                "size": 6,
                "help": "The quantum estimator for the kinetic stress tensor of the physical system.",
                "longhelp": """The quantum estimator for the kinetic stress tensor of the physical system.
                      Returns the 6 independent components in the form [xx, yy, zz, xy, xz, yz].""",
                "func": (
                    lambda: self.tensor2vec(
                        self.kstress_cv() / (self.cell.V * self.beads.nbeads)
                    )
                ),
            },
            "virial_cv": {
                "dimension": "pressure",
                "size": 6,
                "help": "The quantum estimator for the virial stress tensor of the physical system.",
                "longhelp": """The quantum estimator for the virial stress tensor of the physical system.
                      Returns the 6 independent components in the form [xx, yy, zz, xy, xz, yz].""",
                "func": (
                    lambda: self.tensor2vec(
                        self.forces.vir / (self.cell.V * self.beads.nbeads)
                    )
                ),
            },
            "displacedpath": {
                "dimension": "undefined",
                "help": "The displaced path end-to-end distribution estimator",
                "longhelp": """This is the estimator for the end-to-end distribution, that can be used to calculate the
                      particle momentum distribution as described in in L. Lin, J. A. Morrone, R. Car and M. Parrinello,
                      105, 110602 (2010), Phys. Rev. Lett. Takes arguments 'ux', 'uy' and 'uz', which are the components of
                      the path opening vector. Also takes an argument 'atom', which can be either an atom label or index
                      (zero based) to specify which species to find the end-to-end distribution estimator for. If not
                      specified, all atoms are used. Note that one atom is computed at a time, and that each path opening
                      operation costs as much as a PIMD step. Returns the average over the selected atoms of the estimator of
                      exp(-U(u)) for each frame.""",
                "func": self.get_linlin,
            },
            "scaledcoords": {
                "dimension": "undefined",
                "help": "The scaled coordinates estimators that can be used to compute energy and heat capacity",
                "longhelp": """Returns the estimators that are required to evaluate the scaled-coordinates estimators
                       for total energy and heat capacity, as described in T. M. Yamamoto,
                       J. Chem. Phys., 104101, 123 (2005). Returns eps_v and eps_v', as defined in that paper.
                       As the two estimators have a different dimensions, this can only be output in atomic units.
                       Takes one argument, 'fd_delta', which gives the value of the finite difference parameter used -
                       which defaults to """
                + str(-self._DEFAULT_FINDIFF)
                + """. If the value of 'fd_delta' is negative,
                       then its magnitude will be reduced automatically by the code if the finite difference error
                       becomes too large.""",
                "func": self.get_yama_estimators,
                "size": 2,
            },
            "kcv_scaledcoords": {
                "dimension": "undefined",
                "help": "The scaled coordinates estimators that can be used to compute energy and heat capacity",
                "longhelp": """Returns the estimators that are required to evaluate the scaled-coordinates estimators
                       for total energy and heat capacity, as described in T. M. Yamamoto,
                       J. Chem. Phys., 104101, 123 (2005). Returns eps_v and eps_v', as defined in that paper.
                       As the two estimators have a different dimensions, this can only be output in atomic units.
                       Takes one argument, 'fd_delta', which gives the value of the finite difference parameter used -
                       which defaults to """
                + str(-self._DEFAULT_FINDIFF)
                + """. If the value of 'fd_delta' is negative,
                       then its magnitude will be reduced automatically by the code if the finite difference error
                       becomes too large.""",
                "func": self.get_kcv_estimators,
                "size": 2,
            },
            "sc_scaledcoords": {
                "dimension": "undefined",
                "help": "The Suzuki-Chin scaled coordinates estimators that can be used to compute energy and heat capacity",
                "longhelp": """Returns the estimators that are required to evaluate the scaled-coordinates estimators
                       for total energy and heat capacity, as described in T. M. Yamamoto,
                       J. Chem. Phys., 104101, 123 (2005). Returns eps_v and eps_v', as defined in that paper.
                       As the two estimators have a different dimensions, this can only be output in atomic units.
                       Takes one argument, 'fd_delta', which gives the value of the finite difference parameter used -
                       which defaults to """
                + str(-self._DEFAULT_FINDIFF)
                + """. If the value of 'fd_delta' is negative,
                       then its magnitude will be reduced automatically by the code if the finite difference error
                       becomes too large.""",
                "func": self.get_scyama_estimators,
                "size": 2,
            },
            "isotope_scfep": {
                "dimension": "undefined",
                "size": 7,
                "func": self.get_isotope_yama,
                "help": "The scaled-coordinates free energy perturbation scaled mass KE estimator.",
                "longhelp": """Returns the (many) terms needed to compute the scaled-coordinates free energy
                      perturbation scaled mass KE estimator (M. Ceriotti, T. Markland, J. Chem. Phys. 138, 014112 (2013)).
                      Takes two arguments, 'alpha' and 'atom', which give the
                      scaled mass parameter and the atom of interest respectively, and default to '1.0' and ''. The
                      'atom' argument can either be the label of a particular kind of atom, or an index (zero based)
                      of a specific atom. This property computes, for each atom in the selection, an estimator for
                      the kinetic energy it would have had if it had the mass scaled by alpha. The 7 numbers output
                      are the average over the selected atoms of the log of the weights <h>, the average of the
                      squares <h**2>, the average of the un-weighted scaled-coordinates kinetic energies  <T_CV>
                      and of the squares <T_CV**2>, the log sum of the weights LW=ln(sum(e**(-h))), the sum of the
                      re-weighted kinetic energies, stored as a log modulus and sign, LTW=ln(abs(sum(T_CV e**(-h))))
                      STW=sign(sum(T_CV e**(-h))). In practice, the best estimate of the estimator can be computed
                      as [sum_i exp(LTW_i)*STW_i]/[sum_i exp(LW_i)]. The other terms can be used to compute diagnostics
                      for the statistical accuracy of the re-weighting process. Note that evaluating this estimator costs
                      as much as a PIMD step for each atom in the list. The elements that are output have different
                      units, so the output can be only in atomic units.""",
            },
            "isotope_tdfep": {
                "dimension": "undefined",
                "size": 7,
                "func": self.get_isotope_thermo,
                "help": "The thermodynamic free energy perturbation scaled mass KE estimator.",
                "longhelp": """Returns the (many) terms needed to compute the thermodynamic free energy
                      perturbation scaled mass KE estimator (M. Ceriotti, T. Markland, J. Chem. Phys. 138, 014112 (2013)).
                      Takes two arguments, 'alpha' and 'atom', which give the
                      scaled mass parameter and the atom of interest respectively, and default to '1.0' and ''. The
                      'atom' argument can either be the label of a particular kind of atom, or an index (zero based)
                      of a specific atom. This property computes, for each atom in the selection, an estimator for
                      the kinetic energy it would have had if it had the mass scaled by alpha. The 7 numbers output
                      are the average over the selected atoms of the log of the weights <h>, the average of the
                      squares <h**2>, the average of the un-weighted scaled-coordinates kinetic energies  <T_CV>
                      and of the squares <T_CV**2>, the log sum of the weights LW=ln(sum(e**(-h))), the sum of the
                      re-weighted kinetic energies, stored as a log modulus and sign, LTW=ln(abs(sum(T_CV e**(-h))))
                      STW=sign(sum(T_CV e**(-h))). In practice, the best estimate of the estimator can be computed
                      as [sum_i exp(LTW_i)*STW_i]/[sum_i exp(LW_i)]. The other terms can be used to compute diagnostics
                      for the statistical accuracy of the re-weighting process. Evaluating this estimator is inexpensive,
                      but typically the statistical accuracy is worse than with the scaled coordinates estimator.
                      The elements that are output have different
                      units, so the output can be only in atomic units.""",
            },
            "isotope_zetatd": {
                "dimension": "undefined",
                "size": 3,
                "func": self.get_isotope_zetatd,
                "help": "Thermodynamic isotope fractionation direct estimator in the form of ratios of partition functions.",
                "longhelp": """Returns the (many) terms needed to directly compute the relative probablity of
                      isotope substitution in two different systems/phases. Takes two arguments, 'alpha' , which gives the
                      scaled mass parameter and default to '1.0', and 'atom', which is the label or index of a type of atoms.
                      The 3 numbers output are 1) the average over the excess spring energy for an isotope atom substitution <spr>,
                      2) the average of the squares of the excess spring energy <spr**2>, and 3) the average of the exponential
                      of excess spring energy <exp(-beta*spr)>""",
            },
            "isotope_zetasc": {
                "dimension": "undefined",
                "size": 3,
                "func": self.get_isotope_zetasc,
                "help": "Scaled-coordinates isotope fractionation direct estimator in the form of ratios of partition functions.",
                "longhelp": """Returns the (many) terms needed to directly compute the relative probablity of
                      isotope substitution in two different systems/phases. Takes four arguments, 'alpha' , which gives the
                      scaled mass parameter and default to '1.0', and 'atom', which is the label or index of a type of atoms.
                      The 3 numbers output are 1) the average over the excess potential energy for scaled coordinates <sc>,
                      2) the average of the squares of the excess potential energy <sc**2>, and 3) the average of the exponential
                      of excess potential energy <exp(-beta*sc)>""",
            },
            "chin_weight": {
                "dimension": "undefined",
                "size": 3,
                "func": self.get_chin_correction,
                "help": "The weighting factor in Suzuki-Chin 4th-order PI expansion.",
                "longhelp": """The 3 numbers output are 1) the logarithm of the weighting factor -beta_P delta H,
                      2) the square of the logarithm, and 3) the weighting factor""",
            },
            "ti_weight": {
                "dimension": "undefined",
                "size": 3,
                "func": self.get_ti_correction,
                "help": "The weighting factor in Takahashi-Imada 4th-order PI expansion.",
                "longhelp": """The 3 numbers output are 1) the logarithm of the weighting factor -beta_P delta H,
                      2) the square of the logarithm, and 3) the weighting factor""",
            },
            "ti_pot": {
                "dimension": "energy",
                "size": 1,
                "func": self.get_ti_term,
                "help": "The correction potential in Takahashi-Imada 4th-order PI expansion.",
                "longhelp": """The correction potential in Takahashi-Imada 4th-order PI expansion.
                             Takes an argument 'atom', which can be either an atom label or index (zero based)
                             to specify which species to find the correction term for. If not specified,
                             all atoms are used.""",
            },
            "isotope_zetatd_4th": {
                "dimension": "undefined",
                "size": 5,
                "func": self.get_isotope_zetatd_4th,
                "help": "4th order thermodynamic isotope fractionation direct estimator in the form of ratios of partition functions.",
                "longhelp": """Returns the (many) terms needed to compute the thermodynamic
                          fourth-order direct estimator. Takes two arguments, 'alpha' , which gives the
                          scaled mass parameter and default to '1.0', and 'atom', which is the label or
                          index of a type of atoms. The 5 numbers output are 1) the average over the
                          excess spring energy for an isotope atom substitution <spr>, 2) the average
                          of the squares of the excess spring energy <spr**2>, and 3) the average of
                          the exponential of excess spring energy <exp(-beta*spr)>, and 4-5) Suzuki-Chin
                          and Takahashi-Imada 4th-order reweighing term""",
            },
            "isotope_zetasc_4th": {
                "dimension": "undefined",
                "size": 5,
                "func": self.get_isotope_zetasc_4th,
                "help": "4th order scaled-coordinates isotope fractionation direct estimator in the form of ratios of partition functions.",
                "longhelp": """Returns the (many) terms needed to compute the scaled-coordinates
                          fourth-order direct estimator. Takes two arguments, 'alpha' , which gives the scaled
                          mass parameter and default to '1.0', and 'atom', which is the label or index of a type
                          of atoms. The 5 numbers output are 1) the average over the excess potential energy for
                          an isotope atom substitution <sc>, 2) the average of the squares of the excess potential
                          energy <sc**2>, and 3) the average of the exponential of excess potential energy
                          <exp(-beta*sc)>, and 4-5) Suzuki-Chin and Takahashi-Imada 4th-order reweighing term""",
            },
        }

    def bind(self, system):
        """Binds the necessary objects from the system to calculate the
        required properties.

        Args:
           system: The System object to be bound.
        """

        self.ensemble = system.ensemble
        self.motion = system.motion
        self.beads = system.beads
        self.nm = system.nm
        self.cell = system.cell
        self.forces = system.forces
        self.simul = system.simul
        # dummy beads and forcefield objects so that we can use scaled and
        # displaced path estimators without changing the simulation bead
        # coordinates
        self.dbeads = system.beads.copy()
        self.dcell = system.cell.copy()
        self.dforces = system.forces.copy(self.dbeads, self.dcell)
        self.fqref = None
        self._threadlock = (
            system._propertylock
        )  # lock to avoid concurrent access and messing up with dbeads

        # self.properties_init()  # Initialize the properties here so that all
        # +all variables are accessible (for example to set
        # +the size of the hamiltonian_weights).

    def __getitem__(self, key):
        """Retrieves the item given by key.

        Note that if the key contains a string (arg1; arg2; ... )
        then it will pass the appropriate positional arguments to the
        calculation function of the property. Note the brackets and
        the semi-colon separators. If instead we have the syntax
        (arg1=val1;arg2; ... ), then the keyword/value pair (arg1,val1)
        will be added to the keyword argument list. The appropriate key word
        arguments will then be passed to the calculation function instead.

        Similarly, if the key contains a string {unit}, then it will take
        the string 'unit' and use it to define the units that the property
        is output in.

        Args:
           key: A string contained in property_dict.

        Returns:
           The property labelled by the keyword key, along with its unit
           keyword, and the argument lists for the function used to calculate
           the property specified by the keyword key.
        """

        (key, unit, arglist, kwarglist) = getall(key)
        pkey = self.property_dict[key]

        # pkey["func"](*arglist,**kwarglist) gives the value of the property
        # in atomic units. unit_to_user() returns the value in the user
        # specified units.
        with self._threadlock:
            value = pkey["func"](*arglist, **kwarglist)
        if "dimension" in pkey:
            dimension = pkey["dimension"]
        else:
            dimension = ""
        return value, dimension, unit

    def tensor2vec(self, tensor):
        """Takes a 3*3 symmetric tensor and returns it as a 1D array,
        containing the elements [xx, yy, zz, xy, xz, yz].
        """

        return np.array(
            [
                tensor[0, 0],
                tensor[1, 1],
                tensor[2, 2],
                tensor[0, 1],
                tensor[0, 2],
                tensor[1, 2],
            ]
        )

    def get_atom_vec(self, prop_vec, atom="", bead="-1"):
        """Gives a vector for one atom.

        Args:
           prop_vec: An array from which to take the atomic vector from.
           atom: The index of the atom for which the vector will
              be output.
           bead: The index of the replica of the atom for which the
              vector will be output. If less than 0, then the centroid is used.
        """

        if atom == "":
            raise IndexError("Must specify the index for atom_vec property")
        atom = int(atom)
        bead = int(bead)
        if atom >= self.beads.natoms:
            raise IndexError(
                "Cannot output atom_vec property as atom index %d is larger than the number of atoms"
                % atom
            )
        if bead >= self.beads.nbeads:
            raise IndexError(
                "Cannot output atom_vec property as bead index %d is larger than the number of beads"
                % bead
            )

        if bead < 0:
            atom_vec = np.zeros(3)
            for b in range(self.beads.nbeads):
                atom_vec += prop_vec[b, 3 * atom : 3 * (atom + 1)]
            return atom_vec / float(self.beads.nbeads)
        else:
            return prop_vec[bead, 3 * atom : 3 * (atom + 1)]

    def get_temp(self, atom="", bead="", nm=""):
        """Calculates the MD kinetic temperature.

        In case where a specie or a set of indices is selected, and there are constraints,
        the result might be incorrect, as the kinetic energy is not necessarily proportional
        to the temperature. Rather than using a scaling factor based on the number of degrees
        of freedom, we add fake momenta to all the fixed components, so that a meaningful
        result will be obtained for each subset of coordinates regardless of the constraints
        imposed on the system.

        Args:
           atom: If given, specifies the atom to give the temperature
              for. If not, then the simulation temperature.
        """

        if len(self.motion.fixatoms) > 0:
            for i in self.motion.fixatoms:
                pi = np.tile(
                    np.sqrt(
                        self.beads.m[i]
                        * Constants.kb
                        * self.ensemble.temp
                        * self.beads.nbeads
                    ),
                    3,
                )
                self.beads.p[:, 3 * i : 3 * i + 3] += pi

        if self.motion.fixcom:
            # Adds a fake momentum to the centre of mass. This is the easiest way
            # of getting meaningful temperatures for subsets of the system when there
            # are fixed components
            M = np.sum(self.beads.m)
            pcm = np.tile(
                np.sqrt(M * Constants.kb * self.ensemble.temp * self.beads.nbeads), 3
            )
            vcm = np.tile(pcm / M, self.beads.natoms)

            self.beads.p += self.beads.m3 * vcm

            # Avoid double counting
            if len(self.motion.fixatoms) > 0:
                for i in self.motion.fixatoms:
                    self.beads.p[:, 3 * i : 3 * i + 3] -= np.multiply(
                        self.beads.m[i], pcm / M
                    )

        kemd, ncount = self.get_kinmd(atom, bead, nm, return_count=True)

        if self.motion.fixcom:
            # Removes the fake momentum from the centre of mass.
            self.beads.p -= self.beads.m3 * vcm

        if len(self.motion.fixatoms) > 0:
            # re-fixes the fix atoms
            for i in self.motion.fixatoms:
                self.beads.p[:, 3 * i : 3 * i + 3] = 0.0

        return 2.0 * kemd / (Constants.kb * 3.0 * float(ncount) * self.beads.nbeads)

    def get_kincv(self, atom=""):
        """Calculates the quantum centroid virial kinetic energy estimator.

        Args:
           atom: If given, specifies the atom to give the kinetic energy
              for. If not, the system kinetic energy is given.
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output kinetic energy as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        f = dstrip(self.forces.f)
        # subtracts centroid
        q = dstrip(self.beads.q).copy()
        qc = dstrip(self.beads.qc)
        for b in range(self.beads.nbeads):
            q[b] -= qc

        # zeroes components that are not requested
        ncount = 0
        for i in range(self.beads.natoms):
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                q[:, 3 * i : 3 * i + 3] = 0.0
            else:
                ncount += 1

        acv = np.dot(q.flatten(), f.flatten())
        acv *= -0.5 / self.beads.nbeads
        acv += ncount * 1.5 * Constants.kb * self.ensemble.temp
        # ~ acv = 0.0
        # ~ ncount = 0
        # ~
        # ~ for i in range(self.beads.natoms):
        # ~ if (atom != "" and iatom != i and latom != self.beads.names[i]):
        # ~ continue
        # ~
        # ~ kcv = 0.0
        # ~ k = 3*i
        # ~ for b in range(self.beads.nbeads):
        # ~ kcv += q[b,k]* f[b,k] + q[b,k+1]* f[b,k+1] + q[b,k+2]* f[b,k+2]
        # ~ kcv *= -0.5/self.beads.nbeads
        # ~ kcv += 1.5*Constants.kb*self.ensemble.temp
        # ~ acv += kcv
        # ~ ncount += 1

        if ncount == 0:
            warning(
                "Couldn't find an atom which matched the argument of kinetic energy, setting to zero.",
                verbosity.medium,
            )

        return acv

    def get_scpottd(self):
        """Calculates the Suzuki-Chin thermodynamic potential energy estimator."""
        v = 0.0
        pots = dstrip(self.forces.pots)
        potssc = dstrip(self.forces.potssc)
        for k in range(self.beads.nbeads):
            if k % 2 == 0:
                v += 2.0 * pots[k] / 3.0 + 2.0 * (potssc[k] + pots[k] / 3.0)
            else:
                v += 4.0 * pots[k] / 3.0 + 2.0 * (potssc[k] - pots[k] / 3.0)
        return v / (k + 1)

    def get_sckinop(self, atom=""):
        """Calculates the Suzuki-Chin quantum centroid virial kinetic energy estimator.

        Args:
           atom: If given, specifies the atom to give the kinetic energy
              for. If not, the system kinetic energy is given.
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output kinetic energy as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        q = dstrip(self.beads.q)
        qc = dstrip(self.beads.qc)
        f = dstrip(self.forces.f)

        acv = 0.0
        ncount = 0
        for i in range(self.beads.natoms):
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            kcv = 0.0
            k = 3 * i
            for b in range(0, self.beads.nbeads, 2):
                kcv += (
                    (q[b, k] - qc[k]) * f[b, k]
                    + (q[b, k + 1] - qc[k + 1]) * f[b, k + 1]
                    + (q[b, k + 2] - qc[k + 2]) * f[b, k + 2]
                )
            kcv *= -0.5 / self.beads.nbeads * 2.0
            kcv += 1.5 * Constants.kb * self.ensemble.temp
            acv += kcv
            ncount += 1

        if ncount == 0:
            warning(
                "Couldn't find an atom which matched the argument of kinetic energy, setting to zero.",
                verbosity.medium,
            )

        return acv

    def get_fqvirial(self, ref="", units=""):
        if self.fqref is None:
            if ref == "":
                self.fqref = np.zeros(3 * self.beads.natoms)
            else:
                self.fqref = np.loadtxt(ref).flatten() * unit_to_internal(
                    "length", units, 1
                )
                if len(self.fqref) != 3 * self.beads.natoms:
                    raise ValueError(
                        "Atom number mismatch in reference file for virial_fq"
                    )
        fq = 0.0
        for b in range(self.beads.nbeads):
            fq += np.dot(self.forces.f[b], self.beads.q[b] - self.fqref)

        return fq * 0.5 / self.beads.nbeads

    def get_sckintd(self, atom=""):
        """Calculates the Suzuki-Chin thermodynamic quantum centroid virial kinetic energy estimator.

        Args:
           atom: If given, specifies the atom to give the kinetic energy
              for. If not, the system kinetic energy is given.
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output kinetic energy as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        q = dstrip(self.beads.q)
        qc = dstrip(self.beads.qc)
        f = dstrip(self.forces.f)
        fsc = dstrip(self.forces.fsc)

        acv = 0.0
        ncount = 0
        for i in range(self.beads.natoms):
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            kcv = 0.0
            k = 3 * i
            for b in range(self.beads.nbeads):
                kcv += (
                    (q[b, k] - qc[k]) * (f + fsc)[b, k]
                    + (q[b, k + 1] - qc[k + 1]) * (f + fsc)[b, k + 1]
                    + (q[b, k + 2] - qc[k + 2]) * (f + fsc)[b, k + 2]
                )
                if b % 2 == 0:
                    kcv -= (
                        2
                        * (self.forces.alpha / self.forces.omegan2 / 9.0)
                        * (
                            f[b, k] * f[b, k] / self.forces.beads.m3[b, k]
                            + f[b, k + 1] * f[b, k + 1] / self.forces.beads.m3[b, k + 1]
                            + f[b, k + 2] * f[b, k + 2] / self.forces.beads.m3[b, k + 2]
                        )
                    )
                else:
                    kcv -= (
                        2
                        * ((1.0 - self.forces.alpha) / self.forces.omegan2 / 9.0)
                        * (
                            f[b, k] * f[b, k] / self.forces.beads.m3[b, k]
                            + f[b, k + 1] * f[b, k + 1] / self.forces.beads.m3[b, k + 1]
                            + f[b, k + 2] * f[b, k + 2] / self.forces.beads.m3[b, k + 2]
                        )
                    )
            kcv *= -0.5 / self.beads.nbeads
            kcv += 1.5 * Constants.kb * self.ensemble.temp
            acv += kcv
            ncount += 1

        if ncount == 0:
            warning(
                "Couldn't find an atom which matched the argument of kinetic energy, setting to zero.",
                verbosity.medium,
            )

        return acv

    def get_kintd(self, atom=""):
        """Calculates the quantum centroid virial kinetic energy estimator.

        Args:
           atom: If given, specifies the atom to give the kinetic energy
              for. If not, the system kinetic energy is given.
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output kinetic energy as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        q = dstrip(self.beads.q)
        m = dstrip(self.beads.m)
        PkT32 = 1.5 * Constants.kb * self.ensemble.temp * self.beads.nbeads

        atd = 0.0
        ncount = 0
        for i in range(self.beads.natoms):
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            ktd = 0.0
            for b in range(1, self.beads.nbeads):
                for j in range(3 * i, 3 * (i + 1)):
                    ktd += (q[b, j] - q[b - 1, j]) ** 2
            for j in range(3 * i, 3 * (i + 1)):
                ktd += (q[self.beads.nbeads - 1, j] - q[0, j]) ** 2

            ktd *= -0.5 * m[i] * self.nm.omegan2 / self.beads.nbeads
            ktd += PkT32
            atd += ktd
            ncount += 1

        if ncount == 0:
            warning(
                "Couldn't find an atom which matched the argument of kinetic energy, setting to zero.",
                verbosity.medium,
            )

        return atd

    def get_sckinpr(self):
        """Calculates the quantum centroid virial kinetic energy estimator.
        """

        spring = self.beads.vpath * self.nm.omegan2 / self.beads.nbeads
        PkT32 = (
            1.5
            * Constants.kb
            * self.ensemble.temp
            * self.beads.nbeads
            * self.beads.natoms
        )
        pots = dstrip(self.forces.pots)
        potssc = dstrip(self.forces.potssc)
        v = 0.0

        for k in range(self.beads.nbeads):
            if k % 2 == 0:
                v += potssc[k] + pots[k] / 3.0
            else:
                v += potssc[k] - pots[k] / 3.0
        v = v / self.beads.nbeads

        return PkT32 - spring + v

    def get_kinmd(self, atom="", bead="", nm="", return_count=False):
        """Calculates the classical kinetic energy of the simulation (p^2/2m)

        Args:
           atom: If given, specifies the atom to give the kinetic energy
              for. If not, the simulation kinetic energy is given.
           bead: If given, compute the classical KE of a single bead.
           nm: If given, compute the classical KE of a single normal mode.
        """

        if bead != "" and nm != "":
            raise ValueError(
                "Cannot specify both NM and bead for classical kinetic energy estimator"
            )
        if atom != "":
            try:
                # iatom gives the index of the atom to be studied
                iatom = int(atom)
                latom = ""
                if iatom >= self.beads.natoms:
                    raise IndexError(
                        "Cannot output kinetic energy as atom index %d is larger than the number of atoms"
                        % iatom
                    )
            except ValueError:
                # here 'atom' is a label rather than an index which is stored in latom
                iatom = -1
                latom = atom

        ibead = -1
        if bead != "":
            try:
                # iatom gives the index of the atom to be studied
                ibead = int(bead)
                if ibead >= self.beads.nbeads:
                    raise IndexError(
                        "Bead index %d is larger than the number of beads" % ibead
                    )
            except ValueError:
                raise ValueError("Bead index is not a valid integer")

        inm = -1
        if nm != "":
            try:
                # iatom gives the index of the atom to be studied
                inm = int(nm)
                if inm >= self.beads.nbeads:
                    raise IndexError(
                        "Normal mode index %d is larger than the number of beads" % inm
                    )
            except ValueError:
                raise ValueError("Normal mode index is not a valid integer")

        pnm = dstrip(self.nm.pnm)
        dm3 = dstrip(self.nm.dynm3)
        p = dstrip(self.beads.p)
        m3 = dstrip(self.beads.m3)
        kmd = 0.0
        ncount = 0

        if ibead > -1:
            nbeads = 1
            for i in range(self.beads.natoms):
                if atom != "" and iatom != i and latom != self.beads.names[i]:
                    continue
                k = 3 * i
                kmd += (
                    p[ibead, k] ** 2 + p[ibead, k + 1] ** 2 + p[ibead, k + 2] ** 2
                ) / (2.0 * m3[ibead, k])
                ncount += 1
        elif inm > -1:
            nbeads = 1
            for i in range(self.beads.natoms):
                if atom != "" and iatom != i and latom != self.beads.names[i]:
                    continue
                k = 3 * i
                kmd += (
                    pnm[inm, k] ** 2 + pnm[inm, k + 1] ** 2 + pnm[inm, k + 2] ** 2
                ) / (2.0 * dm3[inm, k])
                ncount += 1
        else:
            nbeads = self.beads.nbeads
            ncount = 0
            if atom == "":
                kmd = self.nm.kin
                ncount = self.beads.natoms
            else:
                for i in range(self.beads.natoms):
                    if atom != "" and iatom != i and latom != self.beads.names[i]:
                        continue
                    k = 3 * i
                    for b in range(self.beads.nbeads):
                        kmd += (
                            pnm[b, k] ** 2 + pnm[b, k + 1] ** 2 + pnm[b, k + 2] ** 2
                        ) / (2.0 * dm3[b, k])
                    ncount += 1

        if ncount == 0:
            warning(
                "Couldn't find an atom which matched the argument of kinetic energy, setting to zero.",
                verbosity.medium,
            )

        if return_count:
            return kmd / nbeads, ncount
        else:
            return kmd / nbeads

    def get_ktens(self, atom=""):
        """Calculates the quantum centroid virial kinetic energy
        TENSOR estimator.

        Args:
           atom: The index of the atom for which the kinetic energy tensor
              is to be output, or the index of the type of atoms for which
              it should be output.
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output kinetic tensor as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        tkcv = np.zeros((6), float)
        ncount = 0
        for i in range(self.beads.natoms):
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            tkcv += self.get_kij(str(i), str(i))
            ncount += 1

        if ncount == 0:
            warning(
                "Couldn't find an atom which matched the argument of kinetic tensor, setting to zero.",
                verbosity.medium,
            )

        return tkcv

    def get_vcom(self, latom="", bead=-1):
        """Computes the center of mass velocity for the system or for a specified species"""

        if bead < 0:
            p = dstrip(self.beads.pc)
        else:
            p = dstrip(self.beads[bead])
        pcom = np.zeros(3)
        tm = 0
        for i in range(self.beads.natoms):
            if latom != "" and latom != self.beads.names[i]:
                continue

            tm += self.beads.m[i]
            pcom += p[3 * i : 3 * (i + 1)]
        pcom /= tm
        return pcom

    def get_kij(self, ni="0", nj="0"):
        """Calculates the quantum centroid virial kinetic energy
        TENSOR estimator for two possibly different atom indices.

        Args:
           ni: The index of atom i.
           nj: The index of atom j.

        Returns:
           The contribution to the kinetic energy tensor estimator from
           the interactions between atom i and atom j.
        """

        i = int(ni)
        j = int(nj)
        if i >= self.beads.natoms:
            raise IndexError(
                "Cannot output kinetic_ij as atom index %d is larger than the number of atoms"
                % i
            )
        if j >= self.beads.natoms:
            raise IndexError(
                "Cannot output kinetic_ij as atom index %d is larger than the number of atoms"
                % j
            )
        mi = self.beads.m[i]
        mj = self.beads.m[j]
        ai = 3 * i
        aj = 3 * j

        q = dstrip(self.beads.q)
        qc = dstrip(self.beads.qc)
        f = dstrip(self.forces.f)

        # I implement this for the most general case. In practice T_ij = <p_i p_j>/(2sqrt(m_i m_j))
        kcv = np.zeros((6), float)
        for b in range(self.beads.nbeads):
            kcv[0] += (
                mi * (q[b, ai] - qc[ai]) * f[b, aj]
                + mj * (q[b, aj] - qc[aj]) * f[b, ai]
            )  # Txx
            kcv[1] += (
                mi * (q[b, ai + 1] - qc[ai + 1]) * f[b, aj + 1]
                + mj * (q[b, aj + 1] - qc[aj + 1]) * f[b, ai + 1]
            )  # Tyy
            kcv[2] += (
                mi * (q[b, ai + 2] - qc[ai + 2]) * f[b, aj + 2]
                + mj * (q[b, aj + 2] - qc[aj + 2]) * f[b, ai + 2]
            )  # Tzz
            kcv[3] += (
                mi * (q[b, ai] - qc[ai]) * f[b, aj + 1]
                + mj * (q[b, aj + 1] - qc[aj + 1]) * f[b, ai]
            )  # Txy
            kcv[4] += (
                mi * (q[b, ai] - qc[ai]) * f[b, aj + 2]
                + mj * (q[b, aj + 2] - qc[aj + 2]) * f[b, ai]
            )  # Txz
            kcv[5] += (
                mi * (q[b, ai + 1] - qc[ai + 1]) * f[b, aj + 2]
                + mj * (q[b, aj + 2] - qc[aj + 2]) * f[b, ai + 1]
            )  # Tyz

        kcv *= -0.5 / (self.beads.nbeads * 2 * np.sqrt(mi * mj))
        if i == j:
            kcv[0:3] += 0.5 * Constants.kb * self.ensemble.temp

        return kcv

    def get_rg(self, atom=""):
        """Calculates the radius of gyration of the ring polymers.

        Args:
           atom: If given, specifies the atom to give the gyration radius
              for. If not, the system average gyration radius is given.
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output gyration radius as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        q = dstrip(self.beads.q)
        qc = dstrip(self.beads.qc)
        nat = self.beads.natoms
        nb = self.beads.nbeads
        rg_tot = 0.0
        ncount = 0
        for i in range(nat):
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            rg_at = 0.0
            for j in range(nb):
                dq = q[j, 3 * i : 3 * (i + 1)] - qc[3 * i : 3 * (i + 1)]
                rg_at += np.dot(dq, dq)
            ncount += 1
            rg_tot += np.sqrt(rg_at / float(nb))

        if ncount == 0:
            raise IndexError(
                "Couldn't find an atom which matched the argument of r_gyration"
            )

        return rg_tot / float(ncount)

    def kstress_sctd(self):
        """Calculates the quantum centroid virial kinetic stress tensor
        estimator.

        Note that this is not divided by the volume or the number of beads.

        Returns:
           A 3*3 tensor with all the components of the tensor.
        """

        kst = np.zeros((3, 3), float)
        q = dstrip(self.beads.q)
        qc = dstrip(self.beads.qc)
        pc = dstrip(self.beads.pc)
        m = dstrip(self.beads.m)
        fall = dstrip(self.forces.f + self.forces.fsc)
        na3 = 3 * self.beads.natoms

        for b in range(self.beads.nbeads):
            for i in range(3):
                for j in range(i, 3):
                    kst[i, j] -= np.dot(q[b, i:na3:3] - qc[i:na3:3], fall[b, j:na3:3])

        # return the CV estimator MULTIPLIED BY NBEADS -- again for consistency with the virial, kstress_MD, etc...
        for i in range(3):
            kst[i, i] += self.beads.nbeads * (np.dot(pc[i:na3:3], pc[i:na3:3] / m))

        return kst

    def kstress_cv(self):
        """Calculates the quantum centroid virial kinetic stress tensor
        estimator.

        Note that this is not divided by the volume or the number of beads.

        Returns:
           A 3*3 tensor with all the components of the tensor.
        """

        kst = np.zeros((3, 3), float)
        q = dstrip(self.beads.q)
        qc = dstrip(self.beads.qc)
        pc = dstrip(self.beads.pc)
        m = dstrip(self.beads.m)
        fall = dstrip(self.forces.f)
        na3 = 3 * self.beads.natoms

        for b in range(self.beads.nbeads):
            for i in range(3):
                for j in range(i, 3):
                    kst[i, j] -= np.dot(q[b, i:na3:3] - qc[i:na3:3], fall[b, j:na3:3])

        # return the CV estimator MULTIPLIED BY NBEADS -- again for consistency with the virial, kstress_MD, etc...
        for i in range(3):
            kst[i, i] += self.beads.nbeads * (np.dot(pc[i:na3:3], pc[i:na3:3] / m))

        return kst

    def opening(self, bead):
        """Path opening function, used in linlin momentum distribution
        estimator.

        Args:
           bead: The index of the bead to shift.
        """

        return bead / float(self.beads.nbeads) + 0.5 * (1.0 / self.beads.nbeads - 1)

    def get_linlin(self, ux="0", uy="0", uz="0", atom=""):
        """Calculates the end-to-end distribution for a particular path opening
        vector.

        Args:
           ux: The x-component of the path opening vector.
           uy: The y-component of the path opening vector.
           uz: The z-component of the path opening vector.
           atom: If given, specifies the atom to give the kinetic energy
              for. If not, the simulation kinetic energy is given.
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output linlin estimator as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        beta = 1.0 / (self.ensemble.temp * Constants.kb)

        u = np.array([float(ux), float(uy), float(uz)])
        u_size = np.dot(u, u)
        q = dstrip(self.beads.q)
        nat = self.beads.natoms
        nb = self.beads.nbeads
        nx_tot = 0.0
        ncount = 0
        for i in range(nat):
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            mass = self.beads.m[i]
            self.dbeads.q[:] = q
            for b in range(nb):
                self.dbeads.q[b, 3 * i : 3 * (i + 1)] += self.opening(b) * u
            dV = self.dforces.pot - self.forces.pot

            n0 = np.exp(-mass * u_size / (2.0 * beta * Constants.hbar ** 2))
            nx_tot += n0 * np.exp(-dV * beta / float(self.beads.nbeads))
            ncount += 1

        if ncount == 0:
            raise IndexError(
                "Couldn't find an atom which matched the argument of linlin"
            )

        return nx_tot / float(ncount)

    def get_kcv_estimators(self, fd_delta=-_DEFAULT_FINDIFF):
        """Calculates the op beta derivative of the centroid virial kinetic energy estimator for the Suzuki-Chin propagator.

        Args:
           fd_delta: the relative finite difference in temperature to apply in
           computing finite-difference quantities. If it is negative, will be
           scaled down automatically to avoid discontinuities in the potential.
        """

        eps = abs(float(fd_delta))
        beta = 1.0 / (Constants.kb * self.ensemble.temp)
        beta2 = beta ** 2
        qc = dstrip(self.beads.qc)
        q = dstrip(self.beads.q)

        self.dcell.h = self.cell.h
        self.dbeads.q[::2] = self.beads.q[::2] + eps * (q - qc)[::2]

        vir1 = (
            np.dot(((q - qc)[::2]).flatten(), (self.forces.f[::2]).flatten())
            / self.beads.nbeads
            * 2.0
        )
        vir2 = (
            np.dot(
                ((q - qc)[::2]).flatten(),
                ((self.dforces.f - self.forces.f)[::2]).flatten() / eps,
            )
            / self.beads.nbeads
            * 2.0
        )

        eop = (
            1.5 * self.beads.natoms / beta
            - (0.50 * vir1)
            + np.mean(self.forces.pots[::2])
        )

        r3 = 1.5 * self.beads.natoms / beta2
        r4 = 0.5 / beta * (vir1) * 1.50
        r5 = 0.5 / beta * (vir2) * 0.50

        return np.asarray([eop, r3 + r4 + r5])

    def get_yama_estimators(self, fd_delta=-_DEFAULT_FINDIFF):
        """Calculates the quantum scaled coordinate kinetic energy estimator.

        Uses a finite difference method to calculate the estimators
        needed to calculate the energy and heat capacity of the system, as
        shown in Takeshi M. Yamamoto, Journal of Chemical Physics,
        104101, 123 (2005). Returns both eps_v and eps_v' as defined in
        the above article. Note that heat capacity is calculated as
        beta**2*kboltzmann*(<eps_v**2> - <eps_v>**2 - <eps_v'>), and the
        energy of the system as <eps_v>.

        Args:
           fd_delta: the relative finite difference in temperature to apply in
           computing finite-difference quantities. If it is negative, will be
           scaled down automatically to avoid discontinuities in the potential.
        """

        dbeta = abs(float(fd_delta))
        beta = 1.0 / (Constants.kb * self.ensemble.temp)
        self.dcell.h = self.cell.h
        qc = dstrip(self.beads.qc)
        q = dstrip(self.beads.q)
        v0 = self.forces.pot / self.beads.nbeads
        while True:
            splus = np.sqrt(1.0 + dbeta)
            sminus = np.sqrt(1.0 - dbeta)

            for b in range(self.beads.nbeads):
                self.dbeads[b].q = qc * (1.0 - splus) + splus * q[b, :]
            vplus = self.dforces.pot / self.beads.nbeads

            for b in range(self.beads.nbeads):
                self.dbeads[b].q = qc * (1.0 - sminus) + sminus * q[b, :]
            vminus = self.dforces.pot / self.beads.nbeads

            # print "DISPLACEMENT CHECK YAMA db: %e, d+: %e, d-: %e, dd: %e" %(dbeta, (vplus-v0)*dbeta, (v0-vminus)*dbeta, abs((vplus+vminus-2*v0)/(vplus-vminus)))

            if (
                fd_delta < 0
                and abs((vplus + vminus - 2 * v0) / (vplus - vminus))
                > self._DEFAULT_FDERROR
                and dbeta > self._DEFAULT_MINFID
            ):
                if dbeta > self._DEFAULT_MINFID:
                    dbeta *= 0.5
                    info(
                        "Reducing displacement in scaled coordinates estimator",
                        verbosity.low,
                    )
                    continue
                else:
                    warning(
                        "Could not converge displacement for scaled coordinate estimators",
                        verbosity.low,
                    )
                    eps = 0.0
                    eps_prime = 0.0
                    break
            else:
                eps = ((1.0 + dbeta) * vplus - (1.0 - dbeta) * vminus) / (2 * dbeta)
                eps += 0.5 * (3 * self.beads.natoms) / beta

                eps_prime = (
                    (1.0 + dbeta) * vplus + (1.0 - dbeta) * vminus - 2 * v0
                ) / (dbeta ** 2 * beta)
                eps_prime -= 0.5 * (3 * self.beads.natoms) / beta ** 2

                break

        return np.asarray([eps, eps_prime])

    def get_scyama_estimators(self, fd_delta=-_DEFAULT_FINDIFF):
        """Calculates the quantum scaled coordinate suzuki-chin kinetic energy estimator for the Suzuki-Chin propagator.

        Uses a finite difference method to calculate the estimators
        needed to calculate the energy and heat capacity of the system, as
        shown in Takeshi M. Yamamoto, Journal of Chemical Physics,
        104101, 123 (2005). Returns both eps_v and eps_v' as defined in
        the above article. Note that heat capacity is calculated as
        beta**2*kboltzmann*(<eps_v**2> - <eps_v>**2 - <eps_v'>), and the
        energy of the system as <eps_v>.

        Args:
           fd_delta: the relative finite difference in temperature to apply in
           computing finite-difference quantities. If it is negative, will be
           scaled down automatically to avoid discontinuities in the potential.
        """

        dbeta = abs(float(fd_delta))
        beta = 1.0 / (Constants.kb * self.ensemble.temp)
        self.dforces.omegan2 = self.forces.omegan2
        self.dforces.alpha = self.forces.alpha
        self.dcell.h = self.cell.h

        qc = dstrip(self.beads.qc)
        q = dstrip(self.beads.q)

        v0 = (self.forces.pot + self.forces.potsc) / self.beads.nbeads

        while True:
            splus = np.sqrt(1.0 + dbeta)
            sminus = np.sqrt(1.0 - dbeta)

            for b in range(self.beads.nbeads):
                self.dbeads[b].q = qc * (1.0 - splus) + splus * q[b, :]
            vplus = (self.dforces.pot + self.dforces.potsc) / self.beads.nbeads

            for b in range(self.beads.nbeads):
                self.dbeads[b].q = qc * (1.0 - sminus) + sminus * q[b, :]
            vminus = (self.dforces.pot + self.dforces.potsc) / self.beads.nbeads

            if (
                fd_delta < 0
                and abs((vplus + vminus - 2 * v0) / (vplus - vminus))
                > self._DEFAULT_FDERROR
            ):
                if dbeta > self._DEFAULT_MINFID:
                    dbeta *= 0.5
                    info(
                        "Reducing displacement in scaled coordinates estimator",
                        verbosity.low,
                    )
                    continue
                else:
                    warning(
                        "Could not converge displacement for scaled coordinate estimators",
                        verbosity.low,
                    )
                    eps = 0.0
                    eps_prime = 0.0
                    break
            else:
                eps = ((1.0 + dbeta) * vplus - (1.0 - dbeta) * vminus) / (2 * dbeta)
                eps += 0.5 * (3 * self.beads.natoms) / beta

                eps_prime = (
                    (1.0 + dbeta) * vplus + (1.0 - dbeta) * vminus - 2 * v0
                ) / (dbeta ** 2 * beta)
                eps_prime -= 0.5 * (3 * self.beads.natoms) / beta ** 2

                break

        return np.asarray([eps, eps_prime])

    def get_isotope_yama(self, alpha="1.0", atom=""):
        """Gives the components of the yamamoto scaled-mass KE estimator
        for a given atom index.

        Args:
           alpha: m'/m the mass ratio
           atom: the index of the atom to compute the isotope fractionation
              pair for, or a label

        Returns:
           a tuple from which one can reconstruct all that is needed to
           compute the SMKEE, and its statistical accuracy:
           (sum_deltah, sum_ke, log(sum(weights)), log(sum(weight*ke)),
              sign(sum(weight*ke)) )
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output scaled-mass kinetic energy estimator as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        alpha = float(alpha)

        atcv = 0.0
        atcv2 = 0.0
        alogr = 0.0
        alogr2 = 0.0
        law = 0.0
        lawke = 0.0
        sawke = 1.0
        ni = 0

        # strips dependency control since we are not gonna change the true beads in what follows
        q = dstrip(self.beads.q)
        #        f = dstrip(self.forces.f)
        qc = dstrip(self.beads.qc)

        for i in range(self.beads.natoms):
            # selects only the atoms we care about
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            ni += 1

            # arranges coordinate-scaled beads in a auxiliary beads object
            self.dbeads.q[:] = q[:]
            for b in range(self.beads.nbeads):
                self.dbeads.q[b, 3 * i : 3 * (i + 1)] = qc[
                    3 * i : 3 * (i + 1)
                ] + np.sqrt(1.0 / alpha) * (
                    q[b, 3 * i : 3 * (i + 1)] - qc[3 * i : 3 * (i + 1)]
                )

            tcv = 0.0
            for b in range(self.beads.nbeads):
                tcv += np.dot(
                    (
                        self.dbeads.q[b, 3 * i : 3 * (i + 1)]
                        - self.dbeads.qc[3 * i : 3 * (i + 1)]
                    ),
                    self.dforces.f[b, 3 * i : 3 * (i + 1)],
                )
            tcv *= -0.5 / self.beads.nbeads
            tcv += 1.5 * Constants.kb * self.ensemble.temp

            logr = (self.dforces.pot - self.forces.pot) / (
                Constants.kb * self.ensemble.temp * self.beads.nbeads
            )

            atcv += tcv
            atcv2 += tcv * tcv

            alogr += logr
            alogr2 += logr * logr

            # accumulates log averages in a way which preserves accuracy
            if ni == 1:
                law = -logr
            else:
                (law, drop) = logsumlog((law, 1.0), (-logr, 1.0))

            # here we need to take care of the sign of tcv, which might as well be
            # negative... almost never but...
            if ni == 1:
                lawke = -logr + np.log(abs(tcv))
                sawke = np.sign(tcv)
            else:
                (lawke, sawke) = logsumlog(
                    (lawke, sawke), (-logr + np.log(abs(tcv)), np.sign(tcv))
                )

        if ni == 0:
            raise IndexError(
                "Couldn't find an atom which matched the argument of isotope_y"
            )

        return np.asarray(
            [alogr / ni, alogr2 / ni, atcv / ni, atcv2 / ni, law, lawke, sawke]
        )

    def get_isotope_thermo(self, alpha="1.0", atom=""):
        """Gives the components of the thermodynamic scaled-mass KE
        estimator for a given atom index.

        Args:
           alpha: m'/m the mass ratio
           atom: the index of the atom to compute the isotope fractionation
              pair for, or a label

        Returns:
           a tuple from which one can reconstruct all that is needed to
           compute the SMKEE:
           (sum_deltah, sum_ke, log(sum(weights)), log(sum(weight*ke)),
              sign(sum(weight*ke)) )
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output scaled-mass kinetic energy estimator as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        alpha = float(alpha)

        atcv = 0.0
        alogr = 0.0
        atcv2 = 0.0
        alogr2 = 0.0
        law = 0.0
        lawke = 0.0
        sawke = 1.0
        ni = 0

        # strips dependency control since we are not gonna change the true beads in what follows
        q = dstrip(self.beads.q)
        f = dstrip(self.forces.f)
        qc = dstrip(self.beads.qc)

        for i in range(self.beads.natoms):
            # selects only the atoms we care about
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            ni += 1

            spr = 0.0
            for b in range(1, self.beads.nbeads):
                for j in range(3 * i, 3 * (i + 1)):
                    spr += (q[b, j] - q[b - 1, j]) ** 2
            for j in range(3 * i, 3 * (i + 1)):
                spr += (q[self.beads.nbeads - 1, j] - q[0, j]) ** 2

            spr *= 0.5 * self.beads.m[i] * self.nm.omegan2

            # centroid virial contribution from atom i
            tcv = 0.0
            for b in range(self.beads.nbeads):
                tcv += np.dot(
                    (q[b, 3 * i : 3 * (i + 1)] - qc[3 * i : 3 * (i + 1)]),
                    f[b, 3 * i : 3 * (i + 1)],
                )
            tcv *= -0.5 / self.beads.nbeads
            tcv += 1.5 * Constants.kb * self.ensemble.temp

            logr = (
                (alpha - 1)
                * spr
                / (Constants.kb * self.ensemble.temp * self.beads.nbeads)
            )

            atcv += tcv
            atcv2 += tcv * tcv
            alogr += logr
            alogr2 += logr * logr

            # accumulates log averages in a way which preserves accuracy
            if ni == 1:
                law = -logr
            else:
                (law, drop) = logsumlog((law, 1.0), (-logr, 1.0))

            # here we need to take care of the sign of tcv, which might as well be
            # negative... almost never but...
            if ni == 1:
                lawke = -logr + np.log(abs(tcv))
                sawke = np.sign(tcv)
            else:
                (lawke, sawke) = logsumlog(
                    (lawke, sawke), (-logr + np.log(abs(tcv)), np.sign(tcv))
                )

        if ni == 0:
            raise IndexError(
                "Couldn't find an atom which matched the argument of isotope_y"
            )

        return np.asarray(
            [alogr / ni, alogr2 / ni, atcv / ni, atcv2 / ni, law, lawke, sawke]
        )

    def get_isotope_zetatd(self, alpha="1.0", atom=""):
        """Gives the components  to directly compute the relative probablity of
           isotope substitution in two different systems/phases.

        Args:
           alpha: m'/m the mass ratio
           atom: the label or index of the atom to compute the isotope fractionation pair for

        Returns:
           a tuple from which one can reconstruct all that is needed to
           compute the relative probability of isotope substitution:
           (spraverage, spr2average, sprexpaverage)
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output scaled-mass kinetic energy estimator as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        alpha = float(alpha)

        sprsum = 0.0
        sprexpsum = 0.0
        spr2sum = 0.0
        ni = 0

        # strips dependency control since we are not gonna change the true beads in what follows
        q = dstrip(self.beads.q)
        betaP = 1.0 / (Constants.kb * self.ensemble.temp * self.beads.nbeads)

        for i in range(self.beads.natoms):
            # selects only the atoms we care about
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            ni += 1

            spr = 0.0
            for b in range(1, self.beads.nbeads):
                for j in range(3 * i, 3 * (i + 1)):
                    spr += (q[b, j] - q[b - 1, j]) ** 2
            for j in range(3 * i, 3 * (i + 1)):
                spr += (q[self.beads.nbeads - 1, j] - q[0, j]) ** 2

            # spr = 0.5*(alpha-1)*m_H*omegan2*sum {(q_i+1 - q_i)**2}
            spr *= 0.5 * (alpha - 1.0) * self.beads.m[i] * self.nm.omegan2
            spr2 = spr * spr
            sprexp = np.exp(-betaP * spr)

            sprsum += spr
            spr2sum += spr2
            sprexpsum += sprexp

        if ni == 0:
            raise IndexError(
                "Couldn't find an atom which matched the argument of isotope_zetatd"
            )

        spraverage = sprsum / ni
        spr2average = spr2sum / ni
        sprexpaverage = sprexpsum / ni

        return np.asarray([spraverage, spr2average, sprexpaverage])

    def get_isotope_zetasc(self, alpha="1.0", atom=""):
        """Gives the components  to directly compute the relative probablity of
           isotope substitution in two different systems/phases.

        Args:
           alpha: m'/m the mass ratio
           atom: the label or index of the atom to compute the isotope fractionation pair for

        Returns:
           a tuple from which one can reconstruct all that is needed to
           compute the relative probability of isotope substitution using
           scaled coordinates:
           (yamaaverage, yama2average, yamaexpaverage)
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output scaled-mass kinetic energy estimator as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        alpha = float(alpha)
        scalefactor = 1.0 / np.sqrt(alpha)
        betaP = 1.0 / (Constants.kb * self.ensemble.temp * self.beads.nbeads)

        scsum = 0.0
        scexpsum = 0.0
        sc2sum = 0.0
        ni = 0

        qc = dstrip(self.beads.qc)
        q = dstrip(self.beads.q)
        v0 = self.forces.pot
        self.dbeads.q = q

        for i in range(self.beads.natoms):
            # selects only the atoms we care about
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            ni += 1

            for b in range(self.beads.nbeads):
                for j in range(3 * i, 3 * (i + 1)):
                    self.dbeads.q[b, j] = (
                        qc[j] * (1.0 - scalefactor) + scalefactor * q[b, j]
                    )

            sc = self.dforces.pot - v0
            sc2 = sc * sc
            scexp = np.exp(-betaP * sc)

            scsum += sc
            sc2sum += sc2
            scexpsum += scexp

            self.dbeads.q = q

        if ni == 0:
            raise IndexError(
                "Couldn't find an atom which matched the argument of isotope_zetasc"
            )
        return np.asarray([scsum / ni, sc2sum / ni, scexpsum / ni])

    def get_isotope_zetatd_4th(self, alpha="1.0", atom=""):
        """Gives the components to directly compute the relative probablity of
           isotope substitution in two different systems/phases.
           Includes extra terms needed for Suzuki-Chin and Takahashi-Imada
           4th-order reweighing.

        Args:
           alpha: m'/m the mass ratio
           atom: the label or index of the atom to compute the isotope fractionation pair for

        Returns:
           a tuple that contains terms for the computation of isotope fractionation:
           (spraverage, spr2average, sprexpaverage)
           and re-weighting terms for higher-order correction
            (ti_weight, chin_weight)
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output scaled-mass kinetic energy estimator as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        alpha = float(alpha)

        tdsum = 0.0
        tdexpsum = 0.0
        td2sum = 0.0
        chinexpsum = 0.0
        tiexpsum = 0.0
        ni = 0

        # strips dependency control since we are not gonna change the true beads in what follows
        q = dstrip(self.beads.q)
        f = dstrip(self.forces.f)
        #        m3 = dstrip(self.beads.m3)         # flake8 complains 'unused variable'
        #        pots = self.forces.pots            # -//-
        betaP = 1.0 / (self.beads.nbeads * Constants.kb * self.ensemble.temp)

        for i in range(self.beads.natoms):
            # selects only the atoms we care about
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            ni += 1

            spr = 0.0
            for b in range(1, self.beads.nbeads):
                for j in range(3 * i, 3 * (i + 1)):
                    spr += (q[b, j] - q[b - 1, j]) ** 2
            for j in range(3 * i, 3 * (i + 1)):
                spr += (q[self.beads.nbeads - 1, j] - q[0, j]) ** 2
            spr *= 0.5 * (alpha - 1.0) * self.beads.m[i] * self.nm.omegan2

            # Suzuki-Chin correction
            chin = 0.0
            for b in range(1, self.beads.nbeads, 2):
                for j in range(3 * i, 3 * (i + 1)):
                    chin += f[b, j] ** 2
            chin *= (
                (1.0 / alpha - 1.0)
                * 1.0
                / self.beads.m[i]
                * (4.0 / 3.0)
                * (1.0 / 12.0)
                / self.nm.omegan2
            )

            # Takahashi-Imada correction
            ti = 0.0
            for b in range(self.beads.nbeads):
                for j in range(3 * i, 3 * (i + 1)):
                    ti += f[b, j] ** 2
            ti *= (
                (1.0 / alpha - 1.0)
                * 1.0
                / self.beads.m[i]
                * (1.0 / 24.0)
                / self.nm.omegan2
            )

            td = spr
            td2 = td * td
            tdexp = np.exp(-betaP * td)
            chinexp = np.exp(-betaP * (spr + chin))
            tiexp = np.exp(-betaP * (spr + ti))

            tdsum += td
            td2sum += td2
            tdexpsum += tdexp
            chinexpsum += chinexp
            tiexpsum += tiexp

        if ni == 0:
            raise IndexError(
                "Couldn't find an atom which matched the argument of isotope_zetatd"
            )

        return np.asarray(
            [tdsum / ni, td2sum / ni, tdexpsum / ni, tiexpsum / ni, chinexpsum / ni]
        )

    def get_isotope_zetasc_4th(self, alpha="1.0", atom=""):
        """Gives the components  to directly compute the relative probablity of
           isotope substitution in two different systems/phases.
           Includes extra terms needed for Suzuki-Chin and Takahashi-Imada
           4th-order reweighing.

        Args:
           alpha: m'/m the mass ratio
           atom: the label or index of the atom to compute the isotope fractionation pair for

        Returns:
           a tuple that contains terms for the computation of isotope fractionation:
           (scaverage, sc2average, scexpaverage)
           and re-weighting terms for higher-order correction
           (ti_weight, chin_weight)
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output scaled-mass kinetic energy estimator as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        alpha = float(alpha)
        scalefactor = 1.0 / np.sqrt(alpha)
        betaP = 1.0 / (Constants.kb * self.ensemble.temp * self.beads.nbeads)

        scsum = 0.0
        scexpsum = 0.0
        sc2sum = 0.0
        chinexpsum = 0.0
        tiexpsum = 0.0

        ni = 0

        qc = dstrip(self.beads.qc)
        q = dstrip(self.beads.q)
        f = dstrip(self.forces.f)
        v0 = self.forces.pot
        pots = self.forces.pots

        for i in range(self.beads.natoms):
            # selects only the atoms we care about
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            ni += 1

            self.dbeads.q[:] = q
            # shifts beads positions
            for b in range(self.beads.nbeads):
                for j in range(3 * i, 3 * (i + 1)):
                    self.dbeads.q[b, j] = (
                        qc[j] * (1.0 - scalefactor) + scalefactor * q[b, j]
                    )

            # computes the potential term in the scaled coordinates estimator
            sc = self.dforces.pot - v0

            # this is the extra correction from Suzuki-Chin terms in the hamiltonian.
            # first, the part with |F(q)|^2. this is the scaled-coordinates F with mass m'
            # minus the original coordinates with mass m
            df = dstrip(self.dforces.f)
            dpots = self.dforces.pots

            # Suzuki-Chin correction
            chin = 0.0
            for b in range(1, self.beads.nbeads, 2):
                for j in range(3 * i, 3 * (i + 1)):
                    chin += df[b, j] ** 2 / alpha - f[b, j] ** 2
            chin *= 1.0 / self.beads.m[i] * (4.0 / 3.0) * (1.0 / 12.0) / self.nm.omegan2

            # then, this is the odd/even correction term to the potential.
            # here there is just the mass-scaling that enters, as there is no explicit mass
            for b in range(0, self.beads.nbeads, 2):
                chin += ((-dpots[b] + dpots[b + 1]) - (-pots[b] + pots[b + 1])) / 3.0

            # Takahashi-Imada correction
            ti = 0.0
            for b in range(self.beads.nbeads):
                for j in range(3 * i, 3 * (i + 1)):
                    ti += df[b, j] ** 2 / alpha - f[b, j] ** 2
            ti *= 1.0 / self.beads.m[i] * (1.0 / 24.0) / self.nm.omegan2

            sc2 = sc * sc
            scexp = np.exp(-betaP * sc)
            chinexp = np.exp(-betaP * (sc + chin))
            tiexp = np.exp(-betaP * (sc + ti))

            scsum += sc
            sc2sum += sc2
            scexpsum += scexp
            chinexpsum += chinexp
            tiexpsum += tiexp

        self.dbeads.q[:] = q[:]
        if ni == 0:
            raise IndexError(
                "Couldn't find an atom which matched the argument of isotope_zetasc"
            )

        return np.asarray(
            [scsum / ni, sc2sum / ni, scexpsum / ni, tiexpsum / ni, chinexpsum / ni]
        )

    def get_chin_correction(self):

        f = dstrip(self.forces.f)
        m3 = dstrip(self.beads.m3)
        pots = self.forces.pots
        betaP = 1.0 / (self.beads.nbeads * Constants.kb * self.ensemble.temp)

        chin = 0.0

        for j in range(self.beads.natoms * 3):
            for b in range(1, self.beads.nbeads, 2):  # only loops on odd beads
                chin += (f[b, j] ** 2) / m3[b, j]

        chin *= (4.0 / 3.0) * (1.0 / 12.0) / self.nm.omegan2

        for b in range(0, self.beads.nbeads, 2):
            chin += (-pots[b] + pots[b + 1]) / 3.0

        chin *= -betaP
        chin2 = chin ** 2
        chinexp = np.exp(chin)

        return np.asarray([chin, chin2, chinexp])

    def get_ti_correction(self):

        f = dstrip(self.forces.f)
        m3 = dstrip(self.beads.m3)
        #        pots = self.forces.pots    # flake8 complains "unused variable"
        betaP = 1.0 / (self.beads.nbeads * Constants.kb * self.ensemble.temp)

        ti = 0.0

        for j in range(self.beads.natoms * 3):
            for b in range(self.beads.nbeads):
                ti += (f[b, j] ** 2) / m3[b, j]

        ti *= (1.0 / 24.0) / self.nm.omegan2

        ti *= -betaP
        ti2 = ti ** 2
        tiexp = np.exp(ti)

        return np.asarray([ti, ti2, tiexp])

    def get_ti_term(self, atom=""):
        """Calculates the TI correction potential.

        Args:
           atom: If given, specifies the atom to give the TI correction
              for. If not, the system kinetic energy is given.
        """

        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.beads.natoms:
                raise IndexError(
                    "Cannot output kinetic energy as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        f = dstrip(self.forces.f)
        m3 = dstrip(self.beads.m3)
        #        pots = self.forces.pots                                                # flake8 complains "unused variable"
        #        betaP = 1.0 / (self.beads.nbeads * Constants.kb * self.ensemble.temp)  # -//-

        ti = 0.0

        ncount = 0
        for i in range(self.beads.natoms):
            if atom != "" and iatom != i and latom != self.beads.names[i]:
                continue

            for j in range(3 * i, 3 * (i + 1)):
                for b in range(self.beads.nbeads):
                    ti += (f[b, j] ** 2) / m3[b, j]

            ncount += 1

        ti *= (1.0 / 24.0) / self.nm.omegan2 / self.beads.nbeads
        if ncount == 0:
            warning(
                "Couldn't find an atom which matched the argument of TI potential, setting to zero.",
                verbosity.medium,
            )

        return ti


class Trajectories(dobject):

    """A simple class to take care of output of trajectory data.

    Attributes:
       system: The system object from which the position data will be
          obtained.
       fatom: A dummy beads object used so that individual replica trajectories
          can be output.
       traj_dict: A dictionary containing all the trajectories that can be
          output.
    """

    def __init__(self):
        """Initialises a Trajectories object."""

        self.traj_dict = {
            # Note that here we want to return COPIES of the different arrays, so we make sure to make an operation in order not to return a reference.
            "positions": {
                "dimension": "length",
                "help": "The atomic coordinate trajectories. Will print out one file per bead, unless the bead attribute is set by the user.",
                "func": (lambda: 1.0 * self.system.beads.q),
            },
            "velocities": {
                "dimension": "velocity",
                "help": "The velocity trajectories. Will print out one file per bead, unless the bead attribute is set by the user.",
                "func": (lambda: self.system.beads.p / self.system.beads.m3),
            },
            "momenta": {
                "dimension": "momentum",
                "help": "The momentum trajectories. Will print out one file per bead, unless the bead attribute is set by the user.",
                "func": (lambda: 1.0 * self.system.beads.p),
            },
            "forces": {
                "dimension": "force",
                "help": "The force trajectories. Will print out one file per bead, unless the bead attribute is set by the user.",
                "func": (lambda: 1.0 * self.system.forces.f),
            },
            "forces_sc": {
                "dimension": "force",
                "help": "The Suzuki-Chin component of force trajectories. Will print out one file per bead, unless the bead attribute is set by the user.",
                "func": (
                    lambda: 1.0 * self.system.forces.f + 1.0 * self.system.forces.fsc
                ),
            },
            "x_centroid": {
                "dimension": "length",
                "help": "The centroid coordinates.",
                "func": (lambda: 1.0 * self.system.beads.qc),
            },
            "v_centroid": {
                "dimension": "velocity",
                "help": "The centroid velocity.",
                "func": (lambda: self.system.beads.pc / self.system.beads.m3[0]),
            },
            "x_centroid_even": {
                "dimension": "length",
                "help": "The suzuki-chin centroid coordinates.",
                "func": (
                    lambda: 2
                    * np.sum(self.system.beads.q[::2, :], axis=0)
                    / self.system.beads.nbeads
                ),
            },
            "v_centroid_even": {
                "dimension": "velocity",
                "help": "The suzuki-chin centroid velocity.",
                "func": (
                    lambda: 2
                    * np.sum(
                        (self.system.beads.p / self.system.beads.m3)[::2, :], axis=0
                    )
                    / self.system.beads.nbeads
                ),
            },
            "x_centroid_odd": {
                "dimension": "length",
                "help": "The suzuki-chin centroid coordinates.",
                "func": (
                    lambda: 2
                    * np.sum(self.system.beads.q[1::2, :], axis=0)
                    / self.system.beads.nbeads
                ),
            },
            "v_centroid_odd": {
                "dimension": "velocity",
                "help": "The suzuki-chin centroid velocity.",
                "func": (
                    lambda: 2
                    * np.sum(
                        (self.system.beads.p / self.system.beads.m3)[1::2, :], axis=0
                    )
                    / self.system.beads.nbeads
                ),
            },
            "p_centroid": {
                "dimension": "momentum",
                "help": "The centroid momentum.",
                "func": (lambda: 1.0 * self.system.beads.pc),
            },
            "f_centroid": {
                "dimension": "force",
                "help": "The force acting on the centroid.",
                "func": (
                    lambda: np.sum(self.system.forces.f, 0)
                    / float(self.system.beads.nbeads)
                ),
            },
            "kinetic_cv": {
                "dimension": "energy",
                "help": "The centroid virial quantum kinetic energy estimator for each atom, resolved into Cartesian components [xx, yy, zz]",
                "func": self.get_akcv,
            },
            "kinetic_od": {
                "dimension": "energy",
                "help": "The off diagonal elements of the centroid virial quantum kinetic energy tensor [xy, xz, yz]",
                "func": self.get_akcv_od,
            },
            "r_gyration": {
                "dimension": "length",
                "help": "The radius of gyration of the ring polymer, for each atom and resolved into Cartesian components [xx, yy, zz]",
                "func": self.get_rg,
            },
            "extras": {
                "help": """The additional data returned by the client code, printed verbatim. Will print
                             out one file per bead, unless the bead attribute is set by the user.""",
                "func": (lambda: self.system.forces.extras),
            },
            "isotope_zetatd": {
                "dimension": "undefined",
                "help": """Thermodynamic isotope fractionation direct estimator in the form of ratios of partition functions. Takes two arguments, 'alpha' , which gives the
                      scaled mass parameter and default to '1.0', and 'atom', which is the label or index of a type of atoms. All the atoms but the selected ones
                      will have zero output""",
                "func": self.get_isotope_zetatd,
            },
            "isotope_zetasc": {
                "dimension": "undefined",
                "help": """Scaled-coordinates isotope fractionation direct estimator in the form of ratios of partition functions. Takes two arguments, 'alpha' , which gives the
                      scaled mass parameter and default to '1.0', and 'atom', which is the label or index of a type of atoms. All the atoms but the selected ones
                      will have zero output""",
                "func": self.get_isotope_zetasc,
            },
        }

    def bind(self, system):
        """ Binds to a system object to fetch atomic and force data.

        Args:
           system: The system object that will be managed by this Trajectories.
        """

        self.system = system
        # dummy beads and forcefield objects so that we can use scaled and
        # displaced path estimators without changing the simulation bead
        # coordinates
        self.dbeads = system.beads.copy()
        self.dcell = system.cell.copy()
        self.dforces = self.system.forces.copy(self.dbeads, self.dcell)
        self._threadlock = system._propertylock

        if system.beads.nbeads >= 2:
            self.scdbeads = system.beads.copy(system.beads.nbeads // 2)
            self.scdcell = system.cell.copy()
            self.scdforces = self.system.forces.copy(self.scdbeads, self.scdcell)

    def get_akcv(self):
        """Calculates the contribution to the kinetic energy due to each degree
        of freedom.
        """

        rv = np.zeros(self.system.beads.natoms * 3)
        for b in range(self.system.beads.nbeads):
            rv[:] += (
                self.system.beads.q[b] - self.system.beads.qc
            ) * self.system.forces.f[b]
        rv *= -0.5 / self.system.beads.nbeads
        rv += 0.5 * Constants.kb * self.system.ensemble.temp
        return rv

    def get_akcv_od(self):
        """Calculates the "off-diagonal" contribution to the kinetic energy tensor
        due to each atom.
        """

        rv = np.zeros((self.system.beads.natoms, 3))
        # helper arrays to make it more obvious what we are computing
        dq = np.zeros((self.system.beads.natoms, 3))
        f = np.zeros((self.system.beads.natoms, 3))
        for b in range(self.system.beads.nbeads):
            dq[:] = (self.system.beads.q[b] - self.system.beads.qc).reshape(
                (self.system.beads.natoms, 3)
            )
            f[:] = self.system.forces.f[b].reshape((self.system.beads.natoms, 3))
            rv[:, 0] += dq[:, 0] * f[:, 1] + dq[:, 1] * f[:, 0]
            rv[:, 1] += dq[:, 0] * f[:, 2] + dq[:, 2] * f[:, 0]
            rv[:, 2] += dq[:, 1] * f[:, 2] + dq[:, 2] * f[:, 1]
        rv *= 0.5
        rv *= -0.5 / self.system.beads.nbeads

        return rv.reshape(self.system.beads.natoms * 3)

    def get_rg(self):
        """Calculates the radius of gyration of the ring polymers.

        Computes separately the x, y, z contributions so that the actual
        gyration radius can be recovered as sqrt(rx^2+ry^2+rz^2).
        """

        q = dstrip(self.system.beads.q)
        qc = dstrip(self.system.beads.qc)
        nat = self.system.beads.natoms
        nb = self.system.beads.nbeads
        rg = np.zeros(3 * nat)
        for i in range(nb):
            for j in range(nat):
                dq = q[i, 3 * j : 3 * (j + 1)] - qc[3 * j : 3 * (j + 1)]
                rg[3 * j : 3 * (j + 1)] += dq * dq
        return np.sqrt(rg / float(nb))

    def get_isotope_zetatd(self, alpha="1.0", atom=""):
        """Get the thermodynamic isotope ratio direct estimator for each atom.
        output format:
        column 1: exponent of the direct estimator
        column 2: square of the exponent
        column 3: td estimator

        Args:
           alpha: m'/m the mass ratio
        """
        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.system.beads.natoms:
                raise IndexError(
                    "Cannot output scaled-mass kinetic energy estimator as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        alpha = float(alpha)

        nat = self.system.beads.natoms
        nb = self.system.beads.nbeads
        zetatd = np.zeros((nat, 3))
        # strips dependency control since we are not gonna change the true beads in what follows
        q = dstrip(self.system.beads.q)

        for i in range(nat):
            # selects only the atoms we care about
            if atom != "" and iatom != i and latom != self.system.beads.names[i]:
                continue

            for b in range(1, nb):
                for j in range(3 * i, 3 * (i + 1)):
                    zetatd[i, 0] += (q[b, j] - q[b - 1, j]) ** 2
            for j in range(3 * i, 3 * (i + 1)):
                zetatd[i, 0] += (q[nb - 1, j] - q[0, j]) ** 2

            zetatd[i, 0] *= (
                0.5 * (alpha - 1.0) * self.system.beads.m[i] * self.system.nm.omegan2
            )

        zetatd[:, 1] = np.square(zetatd[:, 0])
        zetatd[:, 2] = np.exp(
            -1.0 / (Constants.kb * self.system.ensemble.temp * nb) * zetatd[:, 0]
        )

        return zetatd.reshape(nat * 3)

    def get_isotope_zetasc(self, alpha="1.0", atom=""):
        """Get the scaled-coordinates isotope ratio direct estimator for each atom.

        output format:
        column 1: exponent of the direct estimator
        column 2: square of the exponent
        column 3: sc estimator

        Args:
           alpha: m'/m the mass ratio
        """
        try:
            # iatom gives the index of the atom to be studied
            iatom = int(atom)
            latom = ""
            if iatom >= self.system.beads.natoms:
                raise IndexError(
                    "Cannot output scaled-mass kinetic energy estimator as atom index %d is larger than the number of atoms"
                    % iatom
                )
        except ValueError:
            # here 'atom' is a label rather than an index which is stored in latom
            iatom = -1
            latom = atom

        alpha = float(alpha)
        scalefactor = 1.0 / np.sqrt(alpha)
        beta = 1.0 / (Constants.kb * self.system.ensemble.temp)

        nat = self.system.beads.natoms
        nb = self.system.beads.nbeads
        zetasc = np.zeros((nat, 3))

        qc = dstrip(self.system.beads.qc)
        q = dstrip(self.system.beads.q)
        v0 = self.system.forces.pot / nb
        self.dbeads.q = q

        for i in range(nat):
            # selects only the atoms we care about
            if atom != "" and iatom != i and latom != self.system.beads.names[i]:
                continue

            for b in range(nb):
                for j in range(3 * i, 3 * (i + 1)):
                    self.dbeads.q[b, j] = (
                        qc[j] * (1.0 - scalefactor) + scalefactor * q[b, j]
                    )
            zetasc[i, 0] = self.dforces.pot / nb - v0

            self.dbeads.q = q

        zetasc[:, 1] = np.square(zetasc[:, 0])
        zetasc[:, 2] = np.exp(-1.0 * beta * zetasc[:, 0])

        return zetasc.reshape(nat * 3)

    def __getitem__(self, key):
        """Retrieves the item given by key.

        Note that if the key contains a string (arg1; arg2; ... )
        then it will pass the appropriate positional arguments to the
        calculation function of the property. Note the brackets and
        the semi-colon separators. If instead we have the syntax
        (arg1=val1;arg2; ... ), then the keyword/value pair (arg1,val1)
        will be added to the keyword argument list. The appropriate key word
        arguments will then be passed to the calculation function instead.

        Similarly, if the key contains a string {unit}, then it will take
        the string 'unit' and use it to define the units that the trajectory
        is output in.

        Args:
           key: A string contained in trajectory_dict.

        Returns:
           The trajectory labelled by the keyword key, along with its unit
           keyword, and the argument lists for the function used to calculate
           the trajectory specified by the keyword key.
        """

        (key, unit, arglist, kwarglist) = getall(key)
        pkey = self.traj_dict[key]

        # pkey["func"](*arglist,**kwarglist) gives the value of the trajectory
        # in atomic units. unit_to_user() returns the value in the user
        # specified units.

        with self._threadlock:
            value = pkey["func"](*arglist, **kwarglist)
        if "dimension" in pkey:
            dimension = pkey["dimension"]
        else:
            dimension = ""
        return value, dimension, unit
