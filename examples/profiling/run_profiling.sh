#!/bin/bash
# usage: run_profile.sh [-p] natoms
# run a mock simulation using an natoms box, -p activates yappi profiling

# Default value for natoms
natoms=8
# Default profiler options
profiler_options=""

# Process command-line arguments
while getopts ":p" opt; do
  case $opt in
    p)
      profiler_flag=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done
shift $((OPTIND-1))

# Set natoms if provided as an argument
if [ $# -gt 0 ]; then
    natoms=$1
fi

# Set profiler options if -p flag was used
if [ "$profiler_flag" = true ]; then
    profiler_options="-p --profiler-clock=wall --profiler-output=test-$natoms"
fi

# Create test.xyz file
echo $natoms > test-$natoms.xyz
echo "# dummy input structure" >> test-$natoms.xyz
for (( i=0; i<natoms; i++ ))
do
    echo "H 0 0 0" >> test-$natoms.xyz
done

# Replace init.xyz with test.xyz in input.xml and save as test.xml
sed 's/init\.xyz/test-'$natoms'\.xyz/g' input.xml > test-$natoms.xml

# Run i-pi and i-pi-driver in the background
time i-pi test-$natoms.xml $profiler_options &> test-$natoms.log &
for i in `seq 1 10`
do
  sleep 1
  i-pi-driver -m gas -u -a gas
done

# Wait for all background processes to finish
wait

