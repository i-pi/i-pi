# Makefile for the driver tests
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

.PHONY: all install modules clean

FLAGS=-g -O3 -Wall
#FLAGS=-g -O0 -Wall -fbacktrace
CFLAGS=$(FLAGS)
FFLAGS=$(FLAGS) -ffree-line-length-none -ffixed-line-length-none -Wno-maybe-uninitialized
#FFLAGS=-g -O0 -fcheck=all -Wall -pedantic -ffree-line-length-none


MODULES=distance.f90 LJ.f90 SG.f90 pes/pswater.f90 pes/LEPS.f90 LJPolymer.f90 pes/efield.f90 pes/eckart.f90 pes/doublewell.f90 pes/doublewell_1D.f90 pes/MB.f90 pes/harmonic_bath.f90
PES=pes/zundel.f pes/morse.f pes/qtip4pf.f pes/utility.f pes/ch52008.f pes/h2o_dip_pol.f

OBJECTS=$(MODULES:%.f90=%.o) $(PES:%.f=%.o)
FC=gfortran
CC=gcc

# `make MPI=1` compiles the ffmpi MPI transport into driver.x using the MPI
# Fortran wrapper. The default build needs no MPI compiler (the MPI code is
# guarded by the __MPI preprocessor macro and stripped out).
ifdef MPI
  FC = mpif90
  CPPMPI = -D__MPI
  MPIOBJ = fmpi.o
endif

# shm_open lives in libc on glibc>=2.34 and on macOS, but in librt on older
# glibc; link librt on Linux (always present, harmless otherwise).
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  LDLIBS += -lrt
endif

all: install

sockets.o: sockets.c
	$(CC) $(CFLAGS) -c -o sockets.o sockets.c

driver.x: $(OBJECTS) sockets.o fsockets.o $(MPIOBJ) driver.o | $(OBJECTS)
	$(FC) $(FFLAGS) -o driver.x $^ $(LDLIBS)

install: driver.x
	ln -fs ../drivers/f90/driver.x ../../bin/i-pi-driver

driver_pure.x: $(OBJECTS) fsockets_pure.o driver.o | $(OBJECTS)
	$(FC) $(FFLAGS) -o driver_pure.x $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

LJ.o: distance.o
LJPolymer.o: distance.o
SG.o: distance.o

# driver.o is preprocessed (-cpp) so the __MPI guards are honoured; fmpi.o (the
# MPI transport module) is only built and linked in an MPI build.
fmpi.o: fmpi.f90
	$(FC) $(FFLAGS) -c fmpi.f90 -o fmpi.o

driver.o: driver.f90 LJ.o LJPolymer.o SG.o pes/pswater.o fsockets.o distance.o $(MPIOBJ)
	$(FC) $(FFLAGS) -cpp $(CPPMPI) -c driver.f90 -o driver.o

clean:
	rm -f *.o *.mod *.x */*.mod */*.o ../../bin/i-pi-driver
