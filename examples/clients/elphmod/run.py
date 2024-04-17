#!/usr/bin/env python3

import elphmod
import ipi._driver.driver

driver = elphmod.md.Driver.load("driver.pickle")

ipi._driver.driver.run_driver(unix=True, address="localhost", driver=driver)
