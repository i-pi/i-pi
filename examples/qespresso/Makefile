# Makefile for the CP2K example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean dia-npt-classical dia_big-npt-classical dia-npt-pimd.4 dia_big-npt-pimd.4 h2o li4
all: dia-npt-classical dia_big-npt-classical dia-npt-pimd.4 dia_big-npt-pimd.4 h2o li4

DRIVER:=pw.x
IPI:=i-pi

-include make.in

define run_driver
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(DRIVER) < ../$2 > pw.out --ipi $3 & cd ..; \
	sleep 3; \
  done;
endef

h2o:
	cd h2o; $(IPI) input.xml & sleep 3; \
	$(call run_driver,4,pw.in,localhost:3141) \
	wait

li4:
	cd li4; $(IPI) input.xml & sleep 3; \
	$(call run_driver,8,pw.in,localhost:31415) \
	wait

dia-npt-classical:
	cd dia-npt-classical; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,pw.in, localhost:3142) \
	wait

dia-npt-pimd.4:
	cd dia-npt-pimd.4; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,pw.in, localhost:3143) \
	wait

dia_big-npt-classical:
	cd dia_big-npt-classical; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,pw.in, localhost:3144) \
	wait

dia_big-npt-pimd.4: 
	cd dia_big-npt-pimd.4; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,pw.in, localhost:3145) \

clean:
	rm -rf */\#*\#* */*simulation* */*.out* */*.log* */*run_* */*RESTART* */EXIT
