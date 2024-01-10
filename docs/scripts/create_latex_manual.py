"""Help script which automatically generates the manual files.

Creates a latex file, corresponding to a section of the manual, for each of
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


if not os.path.exists("input_docs"):
    os.mkdir("input_docs")


ref_includes = open("input_ref_idx.tex", "w")
help(xml=True, prefix="manual")
for opt in objects:
    help(
        latex=True, levels=1, option=opt, prefix=("input_docs/" + opt), standalone=False
    )
    ref_includes.write(f"\include{{input_docs/{opt}}}\n")
for opt in list_objects:
    help_list(option=opt, prefix=("input_docs/" + opt), standalone=False)
