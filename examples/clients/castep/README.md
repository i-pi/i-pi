CASTEP - i-PI example
=====================

*******************************************************************************
  NOTE: The example in this directory is not meant to be physically accurate   
  or reasonable. It is just a quick demonstration to show how to get i-pi and
  CASTEP working together. It uses minimal resources so can be run on a modern
  laptop or workstation etc.
*******************************************************************************

Running the CASTEP example
--------------------------

 * As long as you have CASTEP v22 or later, then it will run without any need to
patch or modify the CASTEP source.

 * CASTEP v22 has built-in support for internet sockets. You can specify the socket port
and hostname via the keywords in the CASTEP param file.

 * A simple launching script has been provided. This will run the example if you type:

```bash
./run_me h2o
```

which will then run 10 steps of NVT PIMD for 1 isolated water molecule in a big box. The
script will launch i-pi server and then create a set of sub-directories, 1 per bead, 
copy the CASTEP input files and then launch the separate CASTEP clients.
