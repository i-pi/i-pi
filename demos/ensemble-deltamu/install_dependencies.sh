#!/bin/bash

conda install -c conda-forge rust python=3.10
conda install --strict-channel-priority -c plumed/label/masterclass-2022 -c conda-forge plumed py-plumed

git clone -b move-rascaline git@github.com:bananenpampe/H2O.git

# pip installs
pip install cmake numpy
pip install --extra-index-url https://download.pytorch.org/whl/cpu torch==2.2.0

export RASCALINE_COMMIT=55ad1f14363812639c865801a6d517d8c5f363f6
export METATENSOR_OPERATIONS_COMMIT=646f1f49c1d9cc18c3fa0d62e2c15e630b8d9d8c
export METATENSOR_COMMIT=2248a3c1feb440a32922ffac98a0592474b48f2d

export PIP_FLAGS="--no-deps --no-build-isolation --check-build-dependencies"

pip install $PIP_FLAGS "metatensor @ git+https://github.com/lab-cosmo/metatensor@$METATENSOR_COMMIT"
pip install $PIP_FLAGS "metatensor-core @ git+https://github.com/lab-cosmo/metatensor@$METATENSOR_COMMIT#subdirectory=python/metatensor-core"
pip install $PIP_FLAGS "metatensor-torch @ git+https://github.com/lab-cosmo/metatensor@$METATENSOR_COMMIT#subdirectory=python/metatensor-torch"

pip install $PIP_FLAGS "metatensor-operations @ git+https://github.com/lab-cosmo/metatensor@$METATENSOR_OPERATIONS_COMMIT#subdirectory=python/metatensor-operations"

pip install $PIP_FLAGS "rascaline @ git+https://github.com/luthaf/rascaline@$RASCALINE_COMMIT"
pip install $PIP_FLAGS "rascaline-torch @ git+https://github.com/luthaf/rascaline@$RASCALINE_COMMIT#subdirectory=python/rascaline-torch"

pip install -r requirements.txt


# modifies the i_pi drivers
cp ./H2O/ipi_driver/lightning.py ../../drivers/py/pes/

export PLUMED_KERNEL=$PWD/plumed-masterclass-2022/lib/libplumedKernel.so
export MLIP_DIR=$PWD/H2O/driver/
