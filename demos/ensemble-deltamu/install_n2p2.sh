#!/bin/bash

echo "Cleaning conda cache"

conda clean -all
conda install python=3.10
conda install conda-forge::lammps
conda uninstall plumed
conda install --strict-channel-priority -c plumed/label/ipi-demo -c conda-forge plumed py-plumed

. ../../env.sh
pip install -e ../../.
