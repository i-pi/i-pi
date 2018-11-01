# Makefile for the Lennard-Jones Neon examples
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

.PHONY: all clean high_density mid_density low_density
all: high_density mid_density low_density

driver:=../../../drivers/driver.x
IPI:=../../../i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h $2 -m lj -o 5.270446,1.1663e-4,13.176115 & \
  done;
endef
-include make.in

high_density:
	cd high_density; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,high_density) \
	wait
	
mid_density:
	cd mid_density; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,mid_density) \
	wait

low_density:
	cd low_density; $(IPI) input.xml & sleep 5; \
	$(call run_driver,1,low_density) \
	wait

clean:
	rm -f */lj-* */RESTART */EXIT
