i-PI driver - Python version
============================

This is a Python version of a simple driver code that computes energy and forces and communicates 
them with i-PI using TCP/IP or Unix sockets. 
Run `i-pi-py_driver -h` to see a list of options and a brief explanation. 

Adding a new PES
----------------

Besides a couple of demonstrative potential energy evaluators, this driver is useful to add
potential-energy surface calculator that facilitate interfacing from any Python-based potential. 

Adding a new PES to the driver is easy: take `pes/dummy.py` as a template and add a new 
file to the `pes/` folder. You can look at some of the other drivers to understand how to
process the command line parameters passed through the `-o` option to the driver.

It is important to set the following global variables in the new file:

```python
__DRIVER_NAME__ = "pes"
__DRIVER_CLASS__ = "PES_class"
```

`__DRIVER_NAME__` is the name that will be used to identify the pes in the `-m` option.
`__DRIVER_CLASS__` is the name you use for the class contained in the new PES file.

