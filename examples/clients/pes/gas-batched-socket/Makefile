# Makefile for the batched-socket ideal-gas (no-op) example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean gas
all: gas

driver:=i-pi-py_driver
IPI:=i-pi


gas:
	$(IPI) input.xml & sleep 5; \
        $(driver) -u -a localhost -m gas & \
	wait

clean:
	rm -f *simulation.* RESTART EXIT
