xtb - i-PI example
========================

Run a simulation of the methane molecule with GFN2-xTB.


A version of the tblite code that is compatible with this example can be 
obtained with 

```bash
mamba install tblite-python -c conda-forge
```

the example can be run as usual

```bash
i-pi input.xml &> log.ipi &
i-pi-py_driver -m xtb -o template=ch4.xyz,method=GFN2-xTB -u -a xtb &> log.xtb
```
