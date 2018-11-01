# Makefile for the para-hydrogen examples
# 
# Copyright (C) 2013, Joshua More and Michele Ceriotti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http.//www.gnu.org/licenses/>.

.PHONY: all clean mht-25K mht-20K mht-17K mht-14K nvt rpmd
all: mht-25K mht-20K mht-17K mht-14K nvt rpmd

driver:=../../../drivers/driver.x
IPI:=../../../i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h $2 -m sg -o 15.0 & \
  done
endef
-include make.in

nvt:
	cd nvt; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,ph2_nvt); \
	wait

mht-25K:
	cd MHT-25K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_25K); \
	wait
	
mht-20K:
	cd MHT-20K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_20K); \
	wait

mht-17K:
	cd MHT-17K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_17K); \
	wait

mht-14K:
	cd MHT-14K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_14K); \
	wait

# this rule creates a run_% subdirectory of RPMD, and executes 16 RPMD runs in it,
# always restarting from the previous run with randomized velocities.
RPMD/run_%/test16.out: 
	cd RPMD; mkdir -p run_$*; cp start-$*.chk run_$*/RESTART; \
	for x in `seq 1 16`; do \
	  sed "s/YYY/"$*"/; s/XXX/$$x/;" input.xml > run_$*/input.xml; \
	  $(IPI) run_$*/input.xml & sleep 25; \
	  $(call run_driver,8,run-$*-$$x); wait; \
	  mv RESTART run_$*/; \
	done  
	

rpmd: RPMD/run_1/test16.out RPMD/run_2/test16.out RPMD/run_3/test16.out RPMD/run_4/test16.out
	make -C RPMD scripts; make -C RPMD postproc
	
	
clean:
	rm -f */ph2* */RESTART */EXIT; \
	make -C RPMD clean; rm -rf RPMD/run_*
