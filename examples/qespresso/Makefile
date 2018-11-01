# Makefile for the Quantum Espresso examples
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

.PHONY: all water lithium clean
all: water lithium zundel

PW:=pw.x
CP:=cp.x
IPI:=../../../i-pi

define run_pw
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(PW) < ../$2 > qe.out & cd ..; \
  done;
endef

define run_cp
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(CP) < ../$2 > qe.out & cd ..; \
  done;
endef


-include make.in

water:
	cd h2o; $(IPI) h2o-ipi.xml & sleep 3; \
	$(call run_pw,4,h2o-pw.in) \
	wait

lithium:
	cd li4; $(IPI) li4-ipi.xml & sleep 3; \
	$(call run_pw,4,li4-pw.in) \
	wait

zundel:
	cd zundel-cp; $(IPI) zundel-ipi.xml & sleep 3; \
	$(call run_cp,4,zundel-cp.in) \
	wait
	
clean:
	rm -rf */*_output.* */*.log */run_* */RESTART */EXIT
