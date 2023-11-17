#!/bin/bash

#Set your path to lammps code.
lmp=lmp

#Run equilibrium trajectory.
i-pi input.xml &> output &
sleep 10
$lmp < in.lmp &> out_lammps &
i-pi-driver -u -h h2o-dipole -m water_dip_pol -o 1 &> out
wait

#Let i-pi close all processes and sockets.
sleep 2

#Run noneqm-traj script. Similar to i-pi but requires two inputs: i-pi input and spectra-related input.
./noneqm-traj.py input.xml -e 0.1 -t 500 &> output &
sleep 10
$lmp < in.lmp &> out_lammps &
i-pi-driver -u -h h2o-dipole -m water_dip_pol -o 0 &> out
wait

./noneqm-response.py input.xml -e 0.1 -t 500
