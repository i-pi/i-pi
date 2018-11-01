"""Deals with testing the driver interface."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import nose
from ipi.interfaces.sockets import Driver, InterfaceSocket
from ipi.interfaces.clients import Client, ClientASE


def test_client():
    """Client: startup without socket."""
    Client(_socket=False)


def test_driver():
    """Driver: startup without socket."""
    Driver(socket=None)


def test_interface():
    """InterfaceSocket: startup."""
    InterfaceSocket()


def test_ASE():
    """Socket client for ASE."""

    try:
        from ase import build
        from ase import units
        from ase.calculators.lj import LennardJones
    except ImportError:
        raise nose.SkipTest

    # create ASE atoms and calculator
    atoms = build.bulk('Ar', cubic=True)
    calculator = LennardJones(epsilon=0.997 * units.kJ / units.mol, sigma=3.4, rc=10.0)
    atoms.set_calculator(calculator)

    # try to get potential energy
    atoms.get_potential_energy()

    # create the socket client
    client = ClientASE(atoms, address='ase', _socket=False)
