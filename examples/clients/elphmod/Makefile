# Makefile for the elphmod example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2024 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all model driver run plot clean
all: plot

IPI = i-pi
MODEL = model_hr.dat model.ifc model.epmatwp model.wigner
DRIVER = driver.pickle driver.xyz
RUN = run.out run.pos_0.xyz run.for_0.xyz
PLOT = plot.pdf

model $(MODEL): model.py
	python3 $<

driver $(DRIVER): driver.py $(MODEL)
	python3 $<

run $(RUN): run.py input.xml $(DRIVER)
	$(IPI) $(word 2,$^) &
	sleep 5
	python3 $<
	wait

plot $(PLOT): plot.py $(RUN)
	python3 $<

clean:
	rm -rf $(MODEL) $(DRIVER) $(RUN) $(PLOT) RESTART EXIT '#'*
