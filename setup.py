"""Script which automatically adds the ipi executable to 
/usr/local/pythonX.Y/site-packages.

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
"""

import os
from setuptools import setup, find_packages 

def read(fname):
   return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
   name='i-PI',
   version='1.0',
   description='A Python interface for ab initio path integral molecular dynamics simulations.',
   long_description=read('README.rst'),
   packages=find_packages(),
   scripts= [ 'i-pi' ],
   libraries = [ ('ipi', {'sources': ['drivers/sockets.c']}) ],
   author= "Michele Ceriotti",
   author_email = "michele.ceriotti@gmail.com",
   classifiers = [ 'Development Status :: 5 - Production/Stable' ],
   license='GPLv3'
)
