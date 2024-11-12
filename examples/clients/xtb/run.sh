#! /bin/bash
echo "Running i-PI"
i-pi input.xml &> log.ipi &
sleep 1
echo "Running driver"
i-pi-py_driver -m xtb -o '{"method": "GFN2-xTB", "numbers": [6,1,1,1,1], "periodic": false}' -u -a xtb &> log.xtb 

