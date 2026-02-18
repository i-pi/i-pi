import numpy as np
import sys

import argparse
from ipi.utils.messages import verbosity, info
from ipi.utils.tools import instanton_compute

""" Reads all the information needed from a i-pi RESTART file and compute the partition functions of the reactant, transition state (TS) or
instanton according to J. Phys. Chem. Lett. 7, 437(2016) (Instanton Rate calculations) or J. Chem. Phys. 134, 054109 (2011) (Tunneling Splitting)


Syntax:    python  Instanton_postproc.py  <checkpoint_file> -c <case> -t  <temperature (K)>  (-n <nbeads(full polymer)>) (-freq <freq_reactant.dat>)

Examples for rate calculation:
           python  Instanton_postproc.py   RESTART  -c  instanton    -t   300
           python  Instanton_postproc.py   RESTART  -c  reactant     -t   300            -n 50
           python  Instanton_postproc.py   RESTART  -c    TS         -t   300

Examples for splitting  calculation (2 steps):
         i)   python  Instanton_postproc.py   RESTART  -c  reactant   -t   10  -n 32 --->this generate the 'eigenvalues_reactant.dat' file
         ii)  python  Instanton_postproc.py   RESTART  -c  instanton  -t   10  -freq eigenvalues_reactant.dat



Type python Instanton_postproc.py -h for more information


Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.
"""

# Y. Litman, 2017.

# You can insert the i-pi path with the following lines.
# Uncomment them and adjust the ipi_path variable

# ipi_path='/home/litman/Yair/Instanton/I-PI-mc/i-pi-mc'

# if not (os.path.exists(ipi_path)):
#    print 'We can not find ipi in %s' %ipi_path
#    print 'Please correct the path'
#    sys.exit()
# sys.path.insert(0, ipi_path)

from ipi.engine.simulation import Simulation
from ipi.utils.units import unit_to_internal, Constants
from ipi.utils.instools import red2comp
from ipi.utils.hesstools import clean_hessian
from ipi.engine.motion.instanton import SpringMapper

# UNITS
K2au = unit_to_internal("temperature", "kelvin", 1.0)
kb = Constants.kb
hbar = Constants.hbar
eV2au = unit_to_internal("energy", "electronvolt", 1.0)
cal2au = unit_to_internal("energy", "cal/mol", 1.0)
cm2au = unit_to_internal("frequency", "hertz", 1.0) * 3e10




if __name__ == "__main__":
    # INPUT
    parser = argparse.ArgumentParser(
        description="""Post-processing routine in order to obtain different quantities from an instanton (or instanton related) calculation. These quantities can be used for the calculation of rates or tunneling splittings in the instanton approximation."""
    )
    parser.add_argument("input", help="Restart file")
    parser.add_argument(
        "-c",
        "--case",
        default=False,
        help="Type of the calculation to analyse. Options: 'instanton', 'reactant' or 'TS'.",
    )
    parser.add_argument(
        "-t", "--temperature", type=float, default=0.0, help="Temperature in K."
    )
    parser.add_argument(
        "-asr",
        "--asr",
        default="poly",
        help="Removes the zero frequency vibrational modes depending on the symmerty of the system",
    )
    parser.add_argument(
        "-e", "--energy_shift", type=float, default=0.0, help="Zero of energy in eV"
    )
    parser.add_argument(
        "-f",
        "--filter",
        default=[],
        help="List of atoms indexes to filter (i.e. eliminate its componentes in the position,mass and hessian arrays. It is 0 based.",
        type=int,
        action="append",
    )
    parser.add_argument(
        "-n",
        "--nbeadsR",
        default=0,
        help="Number of beads (full polymer) to compute the approximate partition function (only reactant case)",
        type=int,
    )
    parser.add_argument(
        "-freq",
        "--freq_reac",
        default=None,
        help="List of frequencies of the minimum. Required for splitting calculation.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        default=False,
        action="store_true",
        help="Avoid the Qvib and Qrot calculation in the instanton case.",
    )

    args = parser.parse_args()
    inputt = args.input
    case = args.case
    temp = args.temperature * K2au
    asr = args.asr
    V00 = args.energy_shift
    filt = args.filter
    nbeadsR = args.nbeadsR
    input_freq = args.freq_reac
    quiet = args.quiet
    Verbosity = verbosity
    Verbosity.level = "quiet"

    instanton_compute(
        inputt, case, temp, asr, V00, filt, nbeadsR, input_freq, quiet, Verbosity
    )
