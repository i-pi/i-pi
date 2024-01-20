#!/bin/bash

################################################################
## ******** EXAMPLE SLURM SUBMISSION SCRIPT FOR I-PI ******** ##
################################################################

## These are the usual SLURM specs, not specifically i-PI related
#SBATCH -t 01:00:00
## ^^^ N.B. the job must last a couple of minutes longer 
## than the <total_time> setting inside input.xml

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH --cpus-per-task=1
## ^^^ N.B. it is important that all jobs are allocated to the same
## node because Unix domain sockets use the filesystem

#SBATCH --mem-per-cpu=1GB
## It is good to fix the memory per process, as on some systems
## otherwise the first driver reserves all the RAM

## Needed since SLURM 22.05 so that srun picks up just one CPU per task
export SRUN_CPUS_PER_TASK=1


## ******* Here starts the actual submission script ******* 

## We assume i-pi (and the driver code) are in the path, otherwise
## you have to set this environment variable
IPI_PATH= 

## Input file
IPI_INPUT=input.xml

## Driver command
IPI_DRIVER="i-pi-driver -a slurm-one-node -m zundel -u -v"

export PATH=$PATH:${IPI_PATH}/bin

## Launches i-PI
srun -n 1 i-pi $IPI_INPUT &> log.ipi &

## Gives a few seconds to allow the server to open the Unix socket
## For *very* complicated simulations you may need to increase a bit
sleep 5; 

## Launches the driver code
for nbead in `seq 1 8`; do
    srun -n 1 $IPI_DRIVER &> log.driver.$nbead &
done

## Waits for all jobs to be finished
wait
