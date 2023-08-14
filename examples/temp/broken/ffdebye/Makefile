# Makefile for the CP2K example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all single_particle replay_with_many_beads
all: single_particle replay_with_many_beads

IPI:=i-pi

single_particle:
	cd single_particle; $(IPI) input.xml;
	wait

replay_with_many_beads:
	cd replay_with_many_beads; $(IPI) input.xml;
	wait

clean:
	rm -rf */*simulation* */*RESTART* */*EXIT*
