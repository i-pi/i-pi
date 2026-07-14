# Makefile for the harmonic oscillator example over an inet (TCP) socket
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean harmonic
all: harmonic

driver:=i-pi-py_driver
IPI:=i-pi

# inet sockets use a host address and port (no -u); the Fortran driver
# i-pi-driver -m harmonic -a localhost -p 33333 -o 1 talks the same protocol.
harmonic:
	$(IPI) input.xml & sleep 5; \
	$(driver) -a localhost -p 33333 -m harmonic -o 1 & \
	wait

clean:
	rm -f *simulation.* RESTART EXIT
