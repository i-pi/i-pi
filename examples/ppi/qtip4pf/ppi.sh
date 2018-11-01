#!/usr/local/bin/bash

source ../../../env.sh && wait;
python ../../../tools/py/energies_ppi.py qtip4pf 298 100 electronvolt > energies.log &
python ../../../tools/py/kinetic_energy_ppi.py qtip4pf 298 100 electronvolt > kinetic_energy.log &
python ../../../tools/py/potential_energy_ppi.py qtip4pf 298 100 electronvolt > potential_energy.log &
python ../../../tools/py/total_energy_ppi.py qtip4pf 298 100 electronvolt > total_energy.log &
python ../../../tools/py/effective_temperatures.py qtip4pf 298 100 > effective_temperatures.log &
python ../../../tools/py/rdf_ppi.py qtip4pf 298 O O 50 2.4 3.3 100 > rdf_OO.log &
python ../../../tools/py/rdf_ppi.py qtip4pf 298 H H 50 1.0 1.9 100 > rdf_HH.log &
python ../../../tools/py/rdf_ppi.py qtip4pf 298 O H 50 0.7 1.3 100 > rdf_OH.log &
