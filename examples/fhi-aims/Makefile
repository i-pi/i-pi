# Makefile for the FHI-AIMS examples
# 
# Copyright (C) 2013, Mariana Rossi and Michele Ceriotti
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

.PHONY: all zundel h2o-npt clean
all: zundel h2o-npt

AIMS:=aims.x
IPI:=../../../i-pi

define run_aims
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cp $2_control.in run_$$i/control.in; cp $2_geometry.in run_$$i/geometry.in; cd run_$$i; \
    $(AIMS) > aims.out & cd ..; \
  done;
endef


-include make.in

zundel:
	cd zundel; $(IPI) zundel-ipi.xml & sleep 3; \
	$(call run_aims,3,zundel-aims) \
	wait

h2o-npt:
	cd h2o-npt; $(IPI) h2o-ipi.xml & sleep 3; \
	$(call run_aims,1,h2o-aims) \
	wait	
	
clean:
	rm -rf */*_output.* */*.log */run_* */RESTART */EXIT
	
