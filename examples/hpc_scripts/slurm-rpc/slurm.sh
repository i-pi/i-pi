#!/bin/bash -l
##--------- SLURM SETTINGS TO BE ADJUSTED -----------------

## do not join stdout and stderr
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err

## name of the job
#SBATCH -J SL-RPC

## execute job from the current working directory
#SBATCH -D ./

#SBATCH --ntasks-per-node=10

## Without this, the first Lammps instance assigns all memory to itself,
## prohibiting running others on the same node.
#SBATCH --mem-per-cpu=1GB

#SBATCH --nodes=1

## run time hh:mm:ss
#SBATCH -t 00:10:00
##---------------------------------------------------------

set -e  # stop on error

# Time in seconds for i-PI should be less than Slurm time
# in order to exit smoothly.
NTIME=500

source  <path-to-i-pi-root>/env.sh
IPI_EXE=<path-to-i-pi-root>/bin/i-pi
LMP_EXE=<path-to-lammps-exe>/lmp_mpi

HOST=`echo $HOSTNAME`
echo "IPI NODE ADDRESS IS $HOST"

echo  {"init",$( date )} >> LIST
if [ -f RESTART ]
then
  grep '<step>'  RESTART >> LIST
  # Save the last restart for the case something goes wrong
  cp RESTART RESTART.save
fi

# Remove i-PI stop signal
if [ -f EXIT ]
then
	rm EXIT
fi

############### i-PI ###############
sed -e "s:<address>.*</address>:<address>$HOST</address>:" input.template-inet.xml > INITIAL.xml

# -u key disables buffering of output in python.
# In case of DFT calculations, its overhead is negligible,
# and it makes debugging a lot easier
python3 -u ${IPI_EXE} INITIAL.xml | tee log.ipi &

# Depending on cluster, i-PI may need a few seconds to start
echo "Waiting 10 seconds..."
sleep 10

############## LAMMPS ##############
# It is recommended to run every client in its own subdirectory
# These subdirectories should be prepared in advance
cd full_sys
#  Plain mpi run
#  cp in.template.lmp in.lmp
#  mpirun -np 2 ${LMP_EXE} -in in.lmp &

#  Slurm run
  sed -e "s/localhost/${HOST}/g" in.template.lmp > in.lmp
  srun -n 4 ${LMP_EXE} -in in.lmp -screen log.lmp -log none &
cd ..
cd beads
#  Plain mpi run
#  cp in.template.lmp in.lmp
#  mpirun -np 2 ${LMP_EXE} -in in.lmp &

#  Slurm run
  sed -e "s/localhost/${HOST}/g" in.template.lmp > in.lmp
  srun -n 1 ${LMP_EXE} -in in.lmp -screen log.instance1.lmp -log none &
  srun -n 1 ${LMP_EXE} -in in.lmp -screen log.instance2.lmp -log none &
  srun -n 1 ${LMP_EXE} -in in.lmp -screen log.instance3.lmp -log none &
  srun -n 1 ${LMP_EXE} -in in.lmp -screen log.instance4.lmp -log none &
cd ..
cd contracted
#  Plain mpi run
#  cp in.template.lmp in.lmp
#  mpirun -np 2 ${LMP_EXE} -in in.lmp &

#  Slurm run
  sed -e "s/localhost/${HOST}/g" in.template.lmp > in.lmp
  srun -n 1 ${LMP_EXE} -in in.lmp -screen log.lmp -log none &
cd ..

wait
sleep 2

########### Save files ###########
grep '<step>'  RESTART >> LIST
echo  {"Final and restart",$( date )} >> LIST
echo '1' >> count
l=`cat count|wc -l`

if [ ! -d IPI_LOGS ]
then
    mkdir IPI_LOGS
fi
cp -p log.ipi IPI_LOGS/log.ipi_$l

########## Resubmit ################
# One can resubmit calculations automatically when running
# long simulations. Use it carefully and responsibly,
# otherwise it can lead to endless loop of resubmission of erroneous script.
if [ ! -f "STOP"   ]
    then
        echo "Resubmit: YES" >> LIST
#        sbatch run-slurm.sh
    else
        echo "Stop file is found." >> LIST
fi
