# Makefile for the yaff example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean mil53_ffsocket mil53_ffyaff
all: mil53_ffsocket mil53_ffyaff

IPI:=i-pi

define run_yaff
  for i in `seq 1 $1`; do \
    python run.py & \
  done;
endef
-include make.in

mil53_ffsocket:
	cd mil53_ffsocket; $(IPI) input.xml & sleep 10; \
	$(call run_yaff,4) \
	wait
	
mil53_ffyaff:
	cd mil53_ffyaff; $(IPI) input.xml; \
	wait

clean:
	rm -f */*.out* */*.log */*.restart* */*RESTART* */EXIT */*pos* */*simulation* */yaffdriver.pyc
