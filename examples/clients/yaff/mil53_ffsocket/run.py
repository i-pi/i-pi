#!/usr/bin/env python

import os

from yaff import *

from yaffdriver import YAFFDriver, MySocket

if __name__ == "__main__":
    logf = open("yaff_%s.log" % os.getpid(), "w")
    log._file = logf
    servername = "MIL53"
    system = System.from_file("init.chk")
    ff = ForceField.generate(
        system,
        "pars.txt",
        rcut=15 * angstrom,
        alpha_scale=3.2,
        gcut_scale=1.5,
        smooth_ei=True,
    )

    socket = MySocket(servername, verbose=True)
    driver = YAFFDriver(socket, ff)
    driver.run()
