# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean mht-25K mht-20K mht-17K mht-14K nvt multi rpmd nst-cubic
all: mht-25K mht-20K mht-17K mht-14K nvt multi rpmd nst-cubic

driver:=i-pi-driver
IPI:=i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h $2 -m sg -o 30.0 & \
  done
endef

nvt:
	cd nvt; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,nvt); \
	wait

multi:
	cd multi; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,multi); \
	wait

mht-25K:
	cd mht-25K; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,mht-25K); \
	wait

mht-20K:
	cd mht-20K; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,mht-20K); \
	wait

mht-17K:
	cd mht-17K; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,mht-17K); \
	wait

mht-14K:
	cd mht-14K; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,mht-14K); \
	wait
	
nst-cubic:
	cd nst-cubic; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,nst-cubic); \
	wait

# this rule creates a run_% subdirectory of rpmd, and executes 16 rpmd runs in it,
# always restarting from the previous run with randomized velocities.
rpmd/run_%/test16.out: 
	cd rpmd; mkdir -p run_$*; cp start-$*.chk run_$*/RESTART; \
	for x in `seq 1 16`; do \
	  sed "s/YYY/"$*"/; s/XXX/$$x/;" input.xml > run_$*/input.xml; \
	  $(IPI) run_$*/input.xml & sleep 25; \
	  $(call run_driver,8,run-$*-$$x); wait; \
	  mv RESTART run_$*/; \
	done  

rpmd: rpmd/run_1/test16.out rpmd/run_2/test16.out rpmd/run_3/test16.out rpmd/run_4/test16.out
	make -C rpmd scripts; make -C rpmd postproc

clean:
	rm -f */*simulation* */*RESTART* */*EXIT*; \
	rm -f para*/*_*ph2* para*/PARA*; \
	make -C rpmd clean; rm -rf rpmd/run_*
