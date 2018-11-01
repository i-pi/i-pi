# Makefile for the FHI-AIMS examples
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all zundel h2o-npt clean
all: zundel h2o-npt

AIMS:=aims.x
IPI:=i-pi

define run_aims
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cp $2_control.in run_$$i/control.in; cp $2_geometry.in run_$$i/geometry.in; cd run_$$i; \
    $(AIMS) > aims.out & cd ..; \
  done;
endef

zundel:
	cd zundel; $(IPI) input.xml & sleep 3; \
	$(call run_aims,3,aims) \
	wait

h2o-npt:
	cd h2o-npt; $(IPI) input.xml & sleep 3; \
	$(call run_aims,1,aims) \
	wait	

clean:
	rm -rf */*simulation* */*.log */run_* */*RESTART* */EXIT
