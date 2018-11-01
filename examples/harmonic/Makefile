# Makefile for the harmonic oscillator example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean harmonic
all: harmonic

driver:=i-pi-driver
IPI:=i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h harmonic -m harm -o 1 & \
  done;
endef
-include make.in

harmonic:
	$(IPI) input.xml & sleep 5; \
	$(call run_driver,4) \
	wait

clean:
	rm -f *simulation.* RESTART EXIT
