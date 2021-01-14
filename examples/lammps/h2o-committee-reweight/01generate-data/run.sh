#!/bin/bash

i-pi input.xml > simulation.log &
sleep 1
lmp_mpi -in in-1.lmp > /dev/null & 
sleep 1
lmp_mpi -in in-2.lmp > /dev/null &

