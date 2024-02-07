# Makefile for the VASP example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean h2o dia
all: h2o dia

VASP:=vasp_std
IPI:=i-pi

define run_vasp
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cp INCAR KPOINTS POTCAR POSCAR run_$$i; cd run_$$i; \
    $(VASP) & cd ..; \
  done;
endef

h2o:
	cd water_1_npt_md_vasp; $(IPI) input.xml & sleep 5; \
	$(call run_vasp,1) \
	wait

dia:
	cd dia_npt_pimd_vasp; $(IPI) input.xml & sleep 5; \
	$(call run_vasp,4) \
	wait

clean:
	rm -rf */*simulation* */*out* */*.log* */*run_* */*RESTART* */EXIT
