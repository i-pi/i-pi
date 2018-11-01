# Makefile for the CP2K example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all single_particle
all: single_particle 

IPI:=i-pi

single_particle:
	cd single_particle; $(IPI) input.xml;
	wait

clean:
	rm -rf */*simulation* */*RESTART* */*EXIT*
