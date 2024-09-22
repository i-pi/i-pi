xtb - i-PI example
========================

Run a simulation of the methane molecule with GFN2-xTB.


A version of the tblite code that is compatible with this example can be 
obtained with 

```bash
pip install tblite
```

the example can be run as usual

```bash
i-pi input.xml &> log.ipi &
i-pi-py_driver -m xtb -o '{"method": "GFN2-xTB", "numbers": [6,1,1,1,1], "periodic": false}' -u -a xtb &> log.xtb &
```
