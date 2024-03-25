#!/bin/bash

export OMP_NUM_THREADS=1

ipi=i-pi
driver=i-pi-py_driver

# ----- SET THE PLUMED AND IPI SOURCE PATHS -----
PLUMED_PATH=/Users/matthiaskellner/Documents/PhD/packages/ipi-feature/plumed-2.8.3
IPI_PATH=/Users/matthiaskellner/Documents/PhD/packages/ipi-feature/i-pi

rm /tmp/ipi*
source  $PLUMED_PATH/sourceme.sh

plumed_path=$(which plumed)
echo "Plumed is located at: $plumed_path"

export PLUMED_KERNEL=$PLUMED_PATH/src/lib/libplumedKernel.dylib
export MLIP_DIR=$PWD/H2O/driver/

sleep_time=2

source $IPI_PATH/env.sh 

if test -f "simulation.restart";
then
        echo "continue last simulation"
        ${ipi} simulation.restart >> log.i-pi &
else
        ${ipi} input_rascaline.xml > log.i-pi & 
fi

echo "# i-PI is running"
echo "# Waiting for ${sleep_time} (s) before executing driver"
sleep ${sleep_time}

${driver} -m lightning -a pinning -u -o ./rascaline_potential/example_5.ckpt,template.xyz &