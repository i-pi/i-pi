#!/bin/bash

conda install -c conda-forge rust python=3.10

#purge conda cache
conda clean --all

#purge the pip cache
pip cache purge

#install plumed and pyplumed from i-pi-demo channel with crystallization enabled
conda install --strict-channel-priority -c plumed/label/ipi-demo -c conda-forge plumed py-plumed

git clone -b i-pi-demo https://github.com/bananenpampe/H2O.git

# pip installs
pip install cmake numpy
pip install --extra-index-url https://download.pytorch.org/whl/cpu torch==2.3.0

pip install -r requirements.txt

pip install metatensor
pip install metatensor-core
pip install metatensor-operations
pip install metatensor-torch

pip install git+https://github.com/Luthaf/rascaline
pip install --extra-index-url https://download.pytorch.org/whl/cpu git+https://github.com/luthaf/rascaline#subdirectory=python/rascaline-torch

# modifies the i_pi drivers
cp ./H2O/ipi_driver/lightning.py ../../drivers/py/pes/

export MLIP_DIR=$PWD/H2O/driver/

. ../../env.sh
pip install -e ../../.
