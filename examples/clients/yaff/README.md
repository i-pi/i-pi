 -- Examples of i-PI working with Yaff --
 
 -- Breathing transition of MIL-53(Al) --

 * NsT simulation of MIL-53(Al) (4 beads) starting in the large pore phase, 
   which collapses at pressures higher than 18 MPa (see e.g. Sven M.J. Rogge,
   Louis Vanduyfhuys, An Ghysels, Michel Waroquier, Toon Verstraelen, Guillaume
   Maurin and Veronique Van Speybroeck, JCTC 11, p. 5583-5597 (2015)) 

 * Uses the force field of Louis Vanduyfhuys, Toon Verstraelen, Matthias
 Vandichel, Michel Waroquier and Veronique Van Speybroeck, JCTC 8, 
 p. 3217-3231 (2012)

 * These simulations are demonstrations, and NOT suited for production runs.
 
 * It is possible to use a Yaff force field without using the ffsocket.
   This is shown in the second example, where the ffyaff class is used.
 
** Runs the examples automatically:

 * First, it is necessary to download and install Yaff. 
 
 Yaff stands for "Yet another force field". It is a pythonic force-field code
 to test-drive new models. The original motivation to develop Yaff was to 
 provide a good reference implementation of the force fields developed at the 
 Center for Molecular Modeling (CMM) at Ghent University.
 
 In its current version, Yaff is general and flexible enough to handle a large
 variety of force field models.
 
 More information, see:
                             https://github.com/molmod/yaff
                             http://molmod.github.io/yaff/ug_install.html

 * Make sure you source the env.sh file in the root directory of i-pi to
   set all the necessary enviroment variables.
 
 * The runs can be done automatically using the Makefile provided. The make 
   targets are self-explanatory. 
   
   To run the first example:

$ make mil53_ffsocket

   To run the second example:

$ make mil53_ffyaff

 * To clean up output files:

$ make clean

** Run the examples manually:

 * Example 1: Go to the mil53_ffsocket directory and run

$ i-pi input.xml
 
   the wrapper will start and sit waiting on the UDS /tmp/ipi.
 
 * Open a separate terminal and run the Yaff code using:
 
$ python run.py

   You can run multiple instances of the code.

 * Example 2: Go to the mil53_ffyaff directory and run

$ i-pi input.xml
