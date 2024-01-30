#!/usr/local/bin/bash

source ../../../env.sh && wait;
python ../../../tools/py/energies_ppi.py lj 10.0 1000 atomic_unit &
python ../../../tools/py/effective_temperatures.py lj 10.0 1000 > effective_temperatures.log &
python ../../../tools/py/rdf_ppi.py lj 10.0 Ne Ne 50 2.2 4 1000 > rdf_log_$folder.log &
