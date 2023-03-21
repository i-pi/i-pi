*******************************************************************************
  NOTE: These examples are provided as a template to use i-pi and VASP, and 
  are designed to run (relatively) quickly on a normal workstation. The 
  parameters of the path integral MD are NOT converged and should not be used 
  verbatim in production calculations. 
*******************************************************************************

                  -- Running the VASP example --
 * In order to run VASP with i-PI, you will need to apply a patch to the VASP source code to add socket support.
   We provide a python3 script `patch_vasp.py` to automate the patching process that is compatible with VASP version >= 5.4:

   python patch_vasp.py vasp.X.Y.Z/src

   where `vasp.X.Y.Z` is the top level directory of the decompressed VASP source code.

 * Alternatively, there are also version-specific patch files in this directory you can use.
   Run the following in the top level directory of VASP, setting the corresponding version number for X.X.X

   patch -p0 < vasp.X.X.X-ipi.patch
   
 * It is then necessary to compile VASP (for other than the available versions please patch it), 
   and create a make.in file containing the path to vasp_std, e.g.

VASP:=~/bin/vasp_std


** Run the examples automatically:
 
 * To run the water example, just type:

$ make h2o

 * To run the lithium example, just type
 
$ make dia

 * To clean up output files:

$ make clean


  
** Run the examples manually:
 
 * In the main example directory run 

$ python ipi-root/bin/i-pi input.xml
 
 * In another terminal, create up to 4 temporary directories to run vasp_std,
 and launch it from within:
 
$ mkdir run_1
$ cp INCAR KPOINTS POSCAR POTCAR run_1
$ cd run_1
$ $VASP
 
 * Repeat several times to get multiple instances
 
 * To run VASP remotely, you should:
    1. make sure that the machine you run i-pi on has a port accessible 
       from outside
    2. change the <interface> option in input.xml so that <address> points 
       to the external network interface (i.e. give a fully-qualified domain 
       name or IP address rather than localhost) and that the <port>
       selected is accessible from outside
    3. change the "IHOST" and "PORT" in INCAR making it point to the same address 
       and port
    4. specify "IBRION=23" in INCAR
    5. copy the VASP inputs to the remote machine (do not forget the pseudos etc.!) 
       and launch n copies of vasp_std remotely
