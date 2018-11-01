"""Help script which automatically generates the manual files.

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

import os
from help import help, objects
from help_list import help_list, list_objects

if not os.path.exists("input_docs"):
   os.mkdir("input_docs")

help(xml=True, prefix="manual")
for opt in objects:
   help(latex=True, levels=1, option=opt, prefix=("input_docs/" + opt), standalone=False)
for opt in list_objects:
   help_list(option=opt, prefix=("input_docs/" + opt), standalone=False)
