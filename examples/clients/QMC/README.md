PET - i-PI example
========================

Run a simulation of the H3 

```bash
i-pi input.xml &> log.ipi &
i-pi-py_driver -m QMC -o parameter_file,ch4.xyz -a pet -u &> log.pet 
```
