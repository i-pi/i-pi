# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean h2o-fixatoms h2o-alchemical h2o-enmfdphonons h2o-fdphonons h2o-geop h2o-multi \
        h2o-piglet.2 h2o-piglet.4 h2o-piglet.8 h2o-pimd.1 h2o-pimd.4 h2o-pimd+mts.4 \
        h2o-pimd+rpc h2o-pimd+sc h2o-pimd+sc+gle.4 h2o-pimd+sc+mts \
        h2o-pimd+sc+rpc h2o-trpmd ice-npt ice-nst isofsc-vapor isofsc-water \
        isof-vapor isof-water NaCl-cellopt paracetamol-geop paracetamol-phonons benz-graph-slc
all: all clean fixatoms h2o-alchemical h2o-enmfdphonons h2o-fdphonons h2o-geop h2o-multi \
     h2o-piglet.2 h2o-piglet.4 h2o-piglet.8 h2o-pimd.1 h2o-pimd.4 h2o-pimd+mts.4 \
     h2o-pimd+rpc h2o-pimd+sc h2o-pimd+sc+gle.4 h2o-pimd+sc+mts \
     h2o-pimd+sc+rpc h2o-trpmd ice-npt ice-nst isofsc-vapor isofsc-water \
     isof-vapor isof-water NaCl-cellopt paracetamol-geop paracetamol-phonons benz-graph-slc

LAMMPS:=lmp_serial
IPI:=i-pi

define run_lammps
  for i in `seq 1 $1`; do \
    $(LAMMPS) < $2 & \
  done;
endef
-include make.in

fixatoms:
	cd fixatoms; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

h2o-alchemical:
	cd h2o-alchemical; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

h2o-enmfdphonons:
	cd h2o-enmfdphonons; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait
	
h2o-fdphonons:
	cd h2o-fdphonons; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait
	
h2o-geop:
	cd h2o-geop; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

h2o-multi:
	cd h2o-multi; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
        wait

h2o-piglet.2:
	cd h2o-piglet.2; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

h2o-piglet.4:
	cd h2o-piglet.4; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

h2o-piglet.8:
	cd h2o-piglet.8; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

h2o-pimd.1:
	cd h2o-pimd.1; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait
	
h2o-pimd.4:
	cd h2o-pimd.4; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

h2o-pimd+mts.4:
	cd h2o-pimd+mts.4; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in_short.lmp) \
        $(call run_lammps,4,in_long.lmp) \
	wait

h2o-pimd+rpc:
	cd h2o-pimd+rpc; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in_short.lmp) \
        $(call run_lammps,4,in_long.lmp) \
	wait
	
h2o-pimd+sc:
	cd h2o-pimd+sc; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait
	
h2o-pimd+sc+gle.4:
	cd h2o-pimd+sc+gle.4; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait
	
h2o-pimd+sc+mts:
	cd h2o-pimd+sc+mts; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

h2o-pimd+sc+rpc:
	cd h2o-pimd+sc+rpc; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in_short.lmp) \
        $(call run_lammps,4,in_long.lmp) \
	wait
	
h2o-trpmd:
	cd h2o-trpmd; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

ice-npt:
	cd ice-npt; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

ice-nst:
	cd ice-nst; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

isofsc-vapor:
	cd isofsc-vapor; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait
	
isofsc-water:
	cd isofsc-water; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

isof-vapor:
	cd isof-vapor; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

isof-water:
	cd isof-water; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

NaCl-cellopt:
	cd NaCl-cellopt; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

paracetamol-geop:
	cd paracetamol-geop; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait
	
paracetamol-phonons:
	cd paracetamol-phonons; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.lmp) \
	wait

benz-graph-slc:
	cd benz-graph-slc; $(IPI) input.xml & sleep 5; \
	$(call run_lammps,1,in.full) \
	$(call run_lammps,1,in.partcontr) \
	$(call run_lammps,8,in.partbeads) \
	wait 

clean:
	rm -f */*RESTART* */EXIT */log.lammps */*simulation*

