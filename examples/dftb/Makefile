# Makefile for the CP2K example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean zundel c18s4
all: zundel c18s4

DFTB:=dftb+
IPI:=i-pi

-include make.in

zundel:
	cd zundel; $(IPI) input.xml & sleep 5; \
	$(DFTB) > dftb.log; \
	wait

c18s4:
	cd c18s4; $(IPI) input.xml & sleep 5; \
	$(DFTB) > dftb.log; \
        wait

clean:
	rm -f */*simulation* */*output* */*.log* */*tmp* */*charges* */*dftb_pin.hsd*  */*RESTART* */EXIT
