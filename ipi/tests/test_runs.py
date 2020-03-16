"""Tests that the Lennard-Jones test case works properly."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from common import SimulationTest


def test_lj_gas():
    ts = SimulationTest(input="../../test/lj/gas/input.xml", driver="../../drivers/driver.x")
    ts.run()
    # Test properties (e.g. latest positions/temperature etc)
