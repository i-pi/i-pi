#! /bin/bash
echo "Running i-PI"
i-pi input.xml &> log.ipi &
sleep 1
echo "Running driver"
i-pi-py_driver -m xtb -o template=ch4.xyz,method=GFN2-xTB -u -a xtb &> log.xtb 

