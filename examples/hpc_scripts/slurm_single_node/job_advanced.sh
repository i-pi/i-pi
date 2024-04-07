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

## Determines the address of the Unix-domain socket used in the input
IPI_ADDRESS=$(grep '<address>' $IPI_INPUT | sed 's/[^>]*>[[:space:]]*//; s/[[:space:]]*<.*//')

## Driver command
IPI_DRIVER="i-pi-driver -a $IPI_ADDRESS -m zundel -u -v"
DRIVER_INPUT=driver.in  # ignored for the example


## We create a UUID to make sure there are no name clashes with other processes,
## or with previous failed runs that left socket files behind
IPI_UUID=$(uuidgen)

## ..and we update the restart file, or create a temporary input that uses this ID
if [ -e RESTART ]; then
    sed -i "s/address>[^<].*</address> $IPI_UUID </" RESTART
else
    sed "s/address>[^<].*</address> $IPI_UUID </" $IPI_INPUT > ${IPI_INPUT}-${IPI_UUID}
fi

## here you may have to update the command-line, or the input driver files
IPI_DRIVER=${IPI_DRIVER//$IPI_ADDRESS/$IPI_UUID}

if [ -e $DRIVER_INPUT ]; then
    ## You may have to adjust the input if the address name is used also 
    ## for something else
    sed "s/$IPI_ADDRESS/$IPI_UUID/" $DRIVER_INPUT > ${DRIVER_INPUT}-${IPI_UUID}
    IPI_DRIVER=${IPI_DRIVER//$DRIVER_INPUT/${DRIVER_INPUT}-${IPI_UUID}}
fi

export PATH=$PATH:${IPI_PATH}/bin

## (Re-)launches i-PI
if [ -e RESTART ]; then
    echo "Restarting i-PI simulation"
    srun -n 1 i-pi RESTART >> log.ipi &
else
    echo "Launching i-PI simulation"
    srun -n 1 i-pi ${IPI_INPUT}-${IPI_UUID} &> log.ipi &
fi

## Gives a few seconds to allow the server to open the Unix socket
## For *very* complicated simulations you may need to increase a bit
sleep 5; 

## Launches the driver code
for nbead in `seq 1 8`; do
    echo "Launching driver ID: $nbead"
    srun -n 1 $IPI_DRIVER &> log.driver.$nbead &
done

## Waits for all jobs to be finished
wait
