# Makefile for the RPMD post-processing files
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

.PHONY: clean scripts postproc

scripts: consolidate vel_est

consolidate: consolidate.f
	gfortran -O3 -o consolidate consolidate.f 
	
vel_est: vel_est.f
	gfortran -O3 -o vel_est vel_est.f

postproc: consolidate vel_est
	for r in `seq 1 4`; do cd run_$$r; for j in `seq 1 16`; do ../vel_est $$j &> /dev/null; done; cd ..; done; \
	echo "4 16" | ./consolidate

clean:
	rm -f consolidate vel_est vel_corr_tot.out */vel_corr*
	
