#!/bin/bash

echo "Cleaning conda cache"

conda clean --all
conda install python=3.10 -y
conda install conda-forge::lammps -y
conda install --no-deps --strict-channel-priority -c plumed/label/ipi-demo -c conda-forge plumed py-plumed -y

. ../../env.sh
pip install -e ../../.