# Makefile for the CP2K example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean npt-classical npt-piglet nvt-classical nvt-piglet nvt-cl
all: npt-classical npt-piglet nvt-classical nvt-piglet nvt-cl

CP2K:=cp2k.x
IPI:=i-pi

define run_cp2k
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(CP2K) -i ../$2 -o cp2k.out & cd ..; \
  done;
endef
-include make.in

npt-classical:
	cd npt-classical; $(IPI) input.xml & sleep 5; \
	$(call run_cp2k,1,in.cp2k) \
	wait

nvt-classical:
	cd nvt-classical; $(IPI) input.xml & sleep 5; \
	$(call run_cp2k,1,in.cp2k) \
	wait

npt-piglet:
	cd npt-piglet; $(IPI) input.xml & sleep 10; \
	$(call run_cp2k,4,in.cp2k) \
	wait

nvt-piglet:
	cd nvt-piglet; $(IPI) input.xml & sleep 10; \
	$(call run_cp2k,4,in.cp2k) \
	wait

#standard output is needed for postprocessing, and saved as "output"
nvt-cl:
	cd nvt-cl; $(IPI) input.xml 2>&1 | tee output & sleep 5; \
	$(call run_cp2k,1,in.cp2k) \
	wait

clean:
	rm -rf */*simulation* */*output* */*.log* */*run_* */*RESTART* */EXIT
