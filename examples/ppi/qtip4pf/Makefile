# Makefile for the qtip4pf example
# 
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.
#
driver:=i-pi-driver
IPI:=i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h qtip4pf -m qtip4pf & \
  done;
endef
-include make.in

all:
	$(IPI) input.xml & sleep 5; \
	$(call run_driver,4) \
	wait

clean:
	rm -f \#* *simulation* *RESTART* EXIT /tmp/ipi_qtip4pf *.log* *.dat*
