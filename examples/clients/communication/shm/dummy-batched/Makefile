# Makefile for the batched-socket ideal-gas (no-op) example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean dummy
all: dummy

driver:=i-pi-driver
IPI:=i-pi


dummy:
	$(IPI) input.xml & sleep 5; \
        $(driver) -u -a localhost -m batch-dummy --shm & \
	wait

clean:
	rm -f *simulation.* RESTART EXIT
