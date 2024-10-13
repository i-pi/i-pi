#! /bin/bash
echo "Running i-PI"
i-pi input.xml &> log.ipi &
sleep 1
echo "Running driver"
<<<<<<< HEAD
i-pi-py_driver -m xtb -o template=ch4.xyz,method=GFN2-xTB -u -a xtb &> log.xtb 
=======
i-pi-py_driver -m xtb -o '{"method": "GFN2-xTB", "numbers": [6,1,1,1,1], "periodic": false}' -u -a xtb &> log.xtb 
>>>>>>> d72ab1f8 (Recovered weird changes wrt main)

