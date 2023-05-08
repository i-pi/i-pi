#!/usr/bin/env python3

import argparse
import sys
from ipi.utils.io import open_backup

try:
    import parmed
    from parmed.periodic_table import Element
except Exception:
    raise RuntimeError(
        """parmed must be installed to use i-pi-amber2xyz. Please run

$~ pip install -U parmed

to install it."""
    )


description = """
Read an Amber parameter and coordinate files and
write a XYZ file in the format expected by i-pi
to initialize a simulation. The Amber coordinate
files can be a formatted restart (.rst7), an
unformatted NetCDF file (.nc or .ncrst), or any
other format supported by the parmed library.
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument("parm", type=str, help="Amber parameter file")
parser.add_argument("crd", type=str, help="Amber coordinate restart file")
parser.add_argument(
    "--filename-out",
    type=str,
    help="Output xyz file. Prints to stdout if not specified.",
)
args = parser.parse_args()


parm = args.parm
crd = args.crd
p = parmed.load_file(parm, xyz=crd)

for a in p.atoms:
    if a.mass < 1:
        raise Exception(f"{parm} atom {a.idx+1} has invalid mass {a.mass}")

if args.filename_out is not None:
    fh = open_backup(args.filename_out, "w")
else:
    fh = sys.stdout

fh.write("%i\n" % (len(p.atoms)))
scell = " ".join(["%12.7f" % (x) for x in p.box])
fh.write("# positions{angstrom} CELL{abcABC}: %s cell{angstrom}\n" % (scell))
for a in p.atoms:
    z = a.atomic_number
    e = Element[z]
    fh.write("%-2s %13.7f %13.7f %13.7f\n" % (e, a.xx, a.xy, a.xz))
