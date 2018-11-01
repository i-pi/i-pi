# Makefile for the harmonic oscillator example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean lj
all: lj

driver:=i-pi-driver
IPI:=i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h lj_4 -m lj -o 5.270446,1.1663e-4,13.176115 & \
  done;
endef
-include make.in

lj:
	source ../../../env.sh; \
	$(IPI) input.xml & sleep 5; \
	$(call run_driver,4) \
	wait

clean:
	rm -f *simulation* *RESTART* EXIT *.log* *.dat*
 
