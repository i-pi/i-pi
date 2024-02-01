"""Help script which automatically generates the manual files.

Creates an rst file, corresponding to a section of the manual, for each of
the classes specified in help.py. It uses help.py to generate information
about the tags for each class, and will include cross-references so that the
title of each tag corresponding to a different class will also be a hyperlink
in the manual to the section corresponding to that class.

Note that any new input class type must be added to the objects
dictionary in help.py and the latex help file must be added to the end of
the manual.lyx file for it to be included in the automatic help generation.

Also creates an xml file with the full list of all the tags.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import os

from help import help, objects
from help_list import help_list, list_objects


if not os.path.exists("input_ref_sections"):
    os.mkdir("input_ref_sections")

ref_includes = open("input_ref_idx.inc", "w")

# help(rst=True, prefix="manual")
for opt in objects:
    help(
        rst=True,
        levels=1,
        option=opt,
        prefix=("input_ref_sections/" + opt),
        standalone=False,
    )
    # rename as .inc to avoid spurious multiple includes
    os.rename(
        "input_ref_sections/" + opt + ".rst", "input_ref_sections/" + opt + ".inc"
    )
    ref_includes.write(f".. include:: input_ref_sections/{opt}.inc\n")

for opt in list_objects:
    help_list(
        option=opt,
        prefix=("input_ref_sections/" + opt),
        latex=False,
        rst=True,
        standalone=False,
    )
    os.rename(
        "input_ref_sections/" + opt + ".rst", "input_ref_sections/" + opt + ".inc"
    )
    ref_includes.write(f".. include:: input_ref_sections/{opt}.inc\n")

# Now I need to reorder some lines
ref_includes.close()

with open("input_ref_idx.inc", "r") as file:
    lines = file.readlines()

# Find and store the relevant lines
property_list_line = next((line for line in lines if "property_list" in line), None)
trajectory_list_line = next((line for line in lines if "trajectory_list" in line), None)

# Remove the original 'property_list' and 'trajectory_list' lines
lines = [
    line for line in lines if line not in [property_list_line, trajectory_list_line]
]

# Insert 'property_list' after 'properties' and 'trajectory_list' after 'trajectory'
for i, line in enumerate(lines):
    if "properties.inc" in line and property_list_line:
        lines.insert(i + 1, property_list_line)
    if "trajectory.inc" in line and trajectory_list_line:
        lines.insert(i + 1, trajectory_list_line)

# Write the modified content back to the file
with open("input_ref_idx.inc", "w") as file:
    file.writelines(lines)
