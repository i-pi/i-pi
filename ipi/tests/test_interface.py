"""Deals with testing the driver interface."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import pytest
try:
    from ase import build
    from ase import units
    from ase.calculators.lj import LennardJones
    has_ase = True
except ImportError:
    has_ase = False


from ipi.interfaces.sockets import Driver, InterfaceSocket
from ipi.interfaces.clients import Client, ClientASE


def test_client():
    """Client: startup without socket."""
    Client(_socket=None)


def test_driver():
    """Driver: startup without socket."""
    Driver(socket=None)


def test_interface():
    """InterfaceSocket: startup."""
    InterfaceSocket()


@pytest.mark.skipif(not has_ase, reason='ASE not installed.')
def test_ASE():
    """Socket client for ASE."""

    # create ASE atoms and calculator
    atoms = build.bulk('Ar', cubic=True)
    calculator = LennardJones(epsilon=0.997 * units.kJ / units.mol, sigma=3.4, rc=10.0)
    atoms.set_calculator(calculator)

    # try to get potential energy
    atoms.get_potential_energy()

    # create the socket client
    ClientASE(atoms, address='ase', _socket=None)
