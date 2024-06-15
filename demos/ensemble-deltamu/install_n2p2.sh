#!/bin/bash

# Check if the miniconda3 directory already exists
if [ ! -d "./miniconda3" ]; then
    echo "Miniconda is not installed. Installing Miniconda..."

    # Create directory for Miniconda
    mkdir -p ./miniconda3

    # Download Miniconda
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./miniconda3/miniconda.sh

    # Make the script executable
    chmod u+rx ./miniconda3/miniconda.sh

    # Execute the installer
    bash ./miniconda3/miniconda.sh -b -u -p ./miniconda3

    # Remove the installer script
    rm -rf ./miniconda3/miniconda.sh

    echo "Miniconda installed successfully."
else
    echo "Miniconda is already installed. Skipping installation."
fi

./miniconda3/bin/conda env remove --name n2p2_install -y
./miniconda3/bin/conda create -n n2p2_install -y

source ./miniconda3/bin/activate n2p2_install

export CONDA_BASE="$PWD/miniconda3/bin"

conda clean --all -y
conda install python=3.10 -y
conda install conda-forge::gsl -y
conda install conda-forge::lammps -y
conda uninstall plumed -y

conda install --no-deps --strict-channel-priority -c plumed/label/ipi-demo -c conda-forge plumed py-plumed -y

. ../../env.sh
pip install -e ../../.
