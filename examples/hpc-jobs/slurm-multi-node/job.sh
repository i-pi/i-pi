#!/bin/bash
################################################################
## ******** EXAMPLE SLURM SUBMISSION SCRIPT FOR I-PI ******** ##
################################################################

## These are the usual SLURM specs, not specifically i-PI related
#SBATCH -t 01:00:00
## ^^^ N.B. the job must last a couple of minutes longer 
## than the <total_time> setting inside input.xml

#SBATCH -N 5
#SBATCH -n 360
#SBATCH --exclusive
## ^^^ N.B. this is a silly setup, in most cases you'd want to run
## highly parallel jobs on each node, while here we have a crazy-cheap
## (and serial) driver. 

##SBATCH --mem-per-cpu=500MB
## It is good to fix the memory per process, as on some systems
## otherwise the first driver reserves all the RAM

## Needed since SLURM 22.05 so that srun picks up the right # of CPUs
export SRUN_CPUS_PER_TASK=1

## ******* Here starts the actual submission script ******* 

## We assume i-pi (and the driver code) are in the path, otherwise
## you have to set this environment variable
IPI_PATH= 

## Input file
IPI_INPUT=input.xml

## Driver command
IPI_DRIVER="i-pi-driver -a HOSTNAME -p PORT -m zundel -v"
## Input file for the driver. Here it is not used, but is still processed to 
## demonstrate the dynamical assignment of the hostname
DRIVER_INPUT=driver.in 

## Determines the address of the Unix-domain socket used in the input
IPI_SRUN=0   # set to 1 if you need to run i-PI on a separate host

if [ -z $IPI_SRUN ]; then
    IPI_ADDRESS=$(hostname)  # fetches the hostname i-PI will run on
else
    IPI_ADDRESS=$(srun -n 1 -N 1 hostname)  # fetches the hostname i-PI will run on
fi

echo "Running on $IPI_ADDRESS"
IPI_PORT=$(grep '<port>' $IPI_INPUT | sed 's/[^>]*>[[:space:]]*//; s/[[:space:]]*<.*//')

## We create a UUID to make sure there are no name clashes with other processes,
## or with previous failed runs that left socket files behind
IPI_UUID=$(uuidgen)

## ..and we update the restart file, or create a temporary input that uses this ID
if [ -e RESTART ]; then
    sed -i "s/address>[^<].*</address> $IPI_ADDRESS </" RESTART
else
    sed "s/address>[^<].*</address> $IPI_ADDRESS </" $IPI_INPUT > ${IPI_INPUT}-${IPI_UUID}
fi

## here you may have to update the command-line, or the input driver files
IPI_DRIVER=${IPI_DRIVER//HOSTNAME/$IPI_ADDRESS}
IPI_DRIVER=${IPI_DRIVER//PORT/$IPI_PORT}

if [ -e $DRIVER_INPUT ]; then
    ## You may have to adjust the input if the address name is used also 
    ## for something else
    sed "s/HOSTNAME/$IPI_ADDRESS/; s/PORT/$IPI_PORT/;" $DRIVER_INPUT > ${DRIVER_INPUT}-${IPI_UUID}
    IPI_DRIVER=${IPI_DRIVER//$DRIVER_INPUT/${DRIVER_INPUT}-${IPI_UUID}}
fi

export PATH=$PATH:${IPI_PATH}/bin

## (Re-)launches i-PI
if [ -e RESTART ]; then
    echo "Restarting i-PI simulation"
    if [ -z $IPI_SRUN ]; then
        i-pi RESTART >> log.ipi &
    else
        srun -n 1 -N 1 i-pi RESTART >> log.ipi &
    fi
else
    echo "Launching i-PI simulation"
    if [ -z $IPI_SRUN ]; then
        i-pi ${IPI_INPUT}-${IPI_UUID} &> log.ipi &
    else
        srun -n 1 -N 1 i-pi ${IPI_INPUT}-${IPI_UUID} &> log.ipi &
    fi
fi

## Gives a few seconds to allow the server to open the Unix socket
## For *very* complicated simulations you may need to increase a bit
sleep 5; 

## Launches the driver code
for nbead in `seq 1 4`; do
    echo "Launching driver ID: $nbead"
    srun -n 1 -N 1 $IPI_DRIVER &> log.driver.$nbead &
done

## Waits for all jobs to be finished
wait
