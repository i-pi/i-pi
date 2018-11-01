#!/usr/local/bin/bash

source ../../../env.sh && wait;
python ../../../tools/py/energies_ppi.py harm 1841.7060385 1000 electronvolt > energies.log &
python ../../../tools/py/effective_temperatures.py harm 1841.7060385 1000 > effective_temperatures.log &
python ../../../tools/py/heat_capacity_ppi.py harm 1841.7060385 1000 > heat_capacity.log &
