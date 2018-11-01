# Makefile for the Lennard-Jones Neon examples
# 
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean high_density mid_density low_density nst
all: high_density mid_density low_density nst

driver:=i-pi-driver
IPI:=i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h $2 -m lj -o 5.270446,1.1663e-4,13.176115 & \
  done;
endef

high_density:
	cd high_density; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,high_density) \
	wait

mid_density:
	cd mid_density; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,mid_density) \
	wait

low_density:
	cd low_density; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,low_density) \
	wait

nst:
	cd nst; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,nst) \
	wait

clean:
	rm -f */*simulation* */*RESTART* */EXIT
