"""Setup script used for package installation."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import os
from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="i-PI",
    version="2.0",
    description="A Python interface for ab initio path integral molecular dynamics simulations.",
    long_description=read("README.rst"),
    packages=find_packages(),
    scripts=["bin/i-pi"],
    author="Michele Ceriotti",
    author_email="michele.ceriotti@gmail.com",
    classifiers=["Development Status :: 5 - Production/Stable"],
    license="GPLv3",
)
