Quantum-ESPRESSO - i-PI example
===============================

*******************************************************************************
  NOTE: These examples are provided as a template to use i-pi and quantum 
  Espresso, and are designed to run (relatively) quickly on a normal 
  workstation. The parameters of both the path integral MD and electronic 
  structure calculation are NOT converged and should not be used verbatim 
  in production calculations.
*******************************************************************************

Running the q-Espresso example
------------------------------

 * First, it is necessary to compile q-e, and create a make.in
   file containing the path to pw.x, e.g.

```bash
PW:=~/bin/pw.x 
```

 * Next, source the env.sh file in the i-pi root.

```bash
source i-pi-root/env.sh
```


** Run the examples automatically:
 
 * To run the water example, just type:

```bash
make h2o
```

 * To run the lithium example, just type

```bash 
make li4
```

 * To run one of the other examples, just type

```bash
make <name of the example>
```

 * To clean up output files:

```bash
make clean
```

  
** Run the examples manually:
 
 * Run i-PI in the main example directory

```bash
python ipi-root/bin/i-pi lithium.xml
```

 * In another terminal, create up to 8 temporary directories to run pw.x,
 and launch it from within:

```bash 
mkdir run_1
cd run_1
$PW --ipi localhost:3141 < ../h2o_pw.in
```

It is important to run client instances in separate folders because QE
uses heavily the disk to store temporary files, which can mess up the 
execution.

 * Repeat several times to get multiple instances
 
* To run q-Espresso remotely, you should:
    1. make sure that the machine you run i-pi on has a port accessible 
       from outside
    2. change the <interface> option in lithium.xml so that <address> points 
       to the external network interface (i.e. give a fully-qualified domain 
       name or IP address rather than localhost) and that the <port>
       selected is accessible from outside
    3. change the `--ipi <address:port> making it point to the same address and port
    4. copy the pw.x inputs to the remote machine (do not forget the pseudos!) 
       and launch n copies of pw.x remotely
