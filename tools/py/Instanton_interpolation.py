"""Instanton_interpolation.py
Reads a hessian file  and/or a positions file (xyz format) and creates an interpolation
that can be used in a further instanton optimization with more beads

Syntax manual:    python  Instanton_interpolation.py -m -xyz <geometry file> -h <hessian file> -n <new-beads(half-polymer)>
Syntax chk:       python  Instanton_interpolation.py -chk  <checkpoint_file>  -n <new-beads(half-polymer)>

Example:   python  Instanton_interpolation.py  -xyz INSTANTON.xyz  -hess  INSTANTON.hess -n 30
           python  Instanton_interpolation.py  -chk RESTART -n 30

Relies on the infrastructure of i-pi, so the ipi package should
be installed in the Python module directory, or the i-pi
main directory must be added to the PYTHONPATH environment variable.

"""

# Y. Litman 2017

import os
import numpy as np
import sys
import argparse

# You can insert the i-pi path with the following lines.
# Uncomment them and adjust the ipi_path variable

# ipi_path='/home/litman/Yair/Instanton/I-PI-mc/i-pi-mc'

# if not (os.path.exists(ipi_path)):
#   print 'We can not find ipi in %s' %ipi_path
#   print 'Please correct the path'
#   sys.exit()
# sys.path.insert(0, ipi_path)

from ipi.utils.io import read_file, print_file
from ipi.utils.nmtransform import nm_rescale
from ipi.utils.units import unit_to_internal
from ipi.utils.tools import interpolate_instanton

if __name__ == "__main__":

    # INPUT
    parser = argparse.ArgumentParser(
        description="""Script for interpolate hessian and/or instanton geometry"""
    )
    parser.add_argument(
        "-m",
        "--manual",
        action="store_true",
        default=False,
        help="Boolean which decides between a checkpoint file or a manual entry.",
    )
    parser.add_argument(
        "-chk",
        "--checkpoint",
        type=str,
        default="None",
        help="Name of the instanton checkpoint file.",
    )
    parser.add_argument(
        "-xyz",
        "--xyz",
        type=str,
        default="None",
        help="Name of the instanton geometry file.",
    )
    parser.add_argument(
        "-hess", "--hessian", type=str, default="None", help="Name of the hessian file."
    )
    parser.add_argument(
        "-n",
        "--nbeadsNew",
        required=True,
        default=0,
        help="New number of beads (half polymer)",
        type=int,
    )

    args = parser.parse_args()
    chk = args.checkpoint
    input_geo = args.xyz
    input_hess = args.hessian
    nbeadsNew = args.nbeadsNew
    manual = args.manual

    if not manual:
        if chk == "None":
            print("Manual mode not specified and checkpoint file name not provided")
            sys.exit()
    else:
        if input_geo == "None":
            print("Manual mode  specified and geometry file name not provided")
            sys.exit()

    interpolate_instanton(chk, input_geo, input_hess, nbeadsNew, manual)
