# Makefile for the harmonic oscillator example
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

.PHONY: all clean harmonic
all: harmonic

driver:=../../drivers/driver.x
IPI:=../../i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h harmonic -m harm -o 1 & \
  done;
endef
-include make.in

harmonic:
	$(IPI) input.xml & sleep 5; \
	$(call run_driver,4) \
	wait

clean:
	rm -f harm.* RESTART EXIT
