# Makefile for the harmonic oscillator example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean harmonic
all: harmonic

driver:=i-pi-py_driver
IPI:=i-pi


harmonic:
	$(IPI) input.xml & sleep 5; \
        $(driver) -u -a localhost -m harmonic -o 1 & \
	wait

clean:
	rm -f *simulation.* RESTART EXIT
