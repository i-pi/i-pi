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

./miniconda3/bin/conda env remove --name rascaline_install -y
./miniconda3/bin/conda create -n rascaline_install -y

source ./miniconda3/bin/activate rascaline_install

export CONDA_BASE="$PWD/miniconda3/bin"

conda install -c conda-forge cxx-compiler gcc rust=1.77 python=3.10 -y

echo $CONDA_PREFIX

export PATH="$CONDA_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

export CC="$CONDA_PREFIX/bin/gcc"
export CXX="$CONDA_PREFIX/bin/g++"

export CMAKE_C_COMPILER="$CONDA_PREFIX/bin/gcc"
export CMAKE_CXX_COMPILER="$CONDA_PREFIX/bin/g++"

#purge conda cache
conda clean --all -y

#purge the pip cache
pip cache purge

#install plumed and pyplumed from i-pi-demo channel with crystallization enabled
conda install --strict-channel-priority -c plumed/label/ipi-demo -c conda-forge plumed py-plumed -y

rm -rf ./H2O/
git clone -b i-pi-demo https://github.com/bananenpampe/H2O.git

# pip installs
pip install cmake==3.29.5.1 numpy==1.26.4
pip install --extra-index-url https://download.pytorch.org/whl/cpu torch==2.3.0

pip install -r requirements.txt

pip install metatensor==0.2.0
pip install metatensor-core==0.1.8
pip install metatensor-operations==0.2.1
pip install metatensor-torch==0.5.1

export RASCALINE_COMMIT=44166158b174ca4162a93238b72dd527347f0800
pip install "rascaline @ git+https://github.com/luthaf/rascaline@$RASCALINE_COMMIT"
pip install --extra-index-url https://download.pytorch.org/whl/cpu "rascaline-torch @ git+https://github.com/luthaf/rascaline@$RASCALINE_COMMIT#subdirectory=python/rascaline-torch" 

# modifies the i_pi drivers
cp ./H2O/ipi_driver/lightning.py ../../drivers/py/pes/

export MLIP_DIR=$PWD/H2O/driver/

. ../../env.sh
pip install -e ../../.
