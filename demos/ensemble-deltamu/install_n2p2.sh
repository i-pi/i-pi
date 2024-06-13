#!/bin/bash

echo "Cleaning conda cache"
export LD_LIBRARY_PATH=

mkdir -p ./miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./miniconda3/miniconda.sh
chmod u+rx ./miniconda3/miniconda.sh
bash ./miniconda3/miniconda.sh -b -u -p ./miniconda3
rm -rf ./miniconda3/miniconda.sh

./miniconda3/bin/conda create -n n2p2_install -y

source ./miniconda3/bin/activate n2p2_install

export CONDA_BASE="$PWD/miniconda3/bin"

conda clean --all
conda install python=3.10 -y
conda install conda-forge::gsl -y
conda install conda-forge::lammps -y
conda uninstall plumed -y

conda install --no-deps --strict-channel-priority -c plumed/label/ipi-demo -c conda-forge plumed py-plumed -y

. ../../env.sh
pip install -e ../../.
