# Makefile for the CP2K example
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

.PHONY: all clean npt-classical npt-piglet nvt-classical nvt-piglet
all: npt-classical npt-piglet nvt-classical nvt-piglet

CP2K:=cp2k.x
IPI:=../../../i-pi

define run_cp2k
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(CP2K) -i ../$2 -o cp2k.out & cd ..; \
  done;
endef
-include make.in

npt-classical:
	cd npt-classical; $(IPI) nptcls_ipi.xml & sleep 5; \
	$(call run_cp2k,1,nptcls_cp2k.in) \
	wait
	
nvt-classical:
	cd nvt-classical; $(IPI) nvtcls_ipi.xml & sleep 5; \
	$(call run_cp2k,1,nvtcls_cp2k.in) \
	wait

npt-piglet:
	cd npt-piglet; $(IPI) nptpgl_ipi.xml & sleep 10; \
	$(call run_cp2k,4,nptpgl_cp2k.in) \
	wait

nvt-piglet:
	cd nvt-piglet; $(IPI) nvtpgl_ipi.xml & sleep 10; \
	$(call run_cp2k,4,nvtpgl_cp2k.in) \
	wait

clean:
	rm -rf */*_output.* */*.log */run_* */RESTART */EXIT
