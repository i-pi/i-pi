#!/usr/bin/env python3

import elphmod

driver = elphmod.md.Driver.load("driver.pickle")

driver.plot(interactive=True, scale=10, pause=0.1)
driver.from_xyz("run.pos_0.xyz")
driver.plot(filename="plot.pdf")
