Examples of usage of the DFTB+ i-PI interface
=============================================

This folder contains examples for running DFTB+ as the backend of i-PI.

Examples exploit fairly advanced features of i-PI, including replica exchange,
threaded execution, and the evaluation of isotope fractionation estimators.

Running the examples
--------------------

Compile DFTB+ and create make.in file containing the path to dftb+ executable, e.g.

```bash
DFTB:=~/bin/dftb+
```

* Run an example automatically, type for instance:

```bash
$ make zundel
```

or

To clean up output files:

```bash
$ make clean
```

* Run an example manually:

In the example directory run

```bash
$ python path/src/i-pi input.xml
```

In another terminal launch dftb+:

```bash
$ path/dftb+ 
```
