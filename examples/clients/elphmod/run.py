#!/usr/bin/env python3

import elphmod
import sys

sys.path.insert(0, "../../../drivers/py")
import driver as ipi_driver

driver = elphmod.md.Driver.load("driver.pickle")

driver.plot(interactive=True)
ipi_driver.run_driver(unix=True, address="localhost", driver=driver)
driver.plot(interactive=False)
