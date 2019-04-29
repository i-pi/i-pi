"""Deals with testing specialisations of HolonomicConstraint."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2019 i-PI developers
# See the "licenses" directory for full license information.


from common import local
from ipi.utils.depend import dstrip
from ipi.engine.initializer import init_file
from ipi.engine.beads import Beads
from ipi.engine.constraints import Replicas, BondLength, BondAngle
import numpy as np
import pytest

def create_beads(fin):
    """Read atom coordinates from file @fin and return an initialised
       instance of beads.
    """
    mode = "xyz"
    ret = init_file(mode, fin)
    atoms = ret[0]
    beads = Beads(atoms[0].natoms, len(atoms))
    for i in range(len(atoms)):
        beads[i] = atoms[i]
    replica_list = []
    for i in range(3*beads.natoms):
        replica_list.append(Replicas(len(atoms), beads, i))

    return beads, replica_list

def bondlength(beads, indices, target):
    """Calculate the mean of the bondlenght over the replicas and
       return the value of the constraint function and its gradient.
    """
    grad = np.empty((beads.nbeads,len(indices)),float)
    for i,idx in enumerate(indices):
        grad[:,i] = dstrip(beads.q[:,idx])
    grad[:,3:] -= grad[:,:3]
    rvec = np.sqrt(np.sum(grad[:,3:]**2, axis=-1))
    sigma = rvec.mean()
    sigma -= target
    grad[:,3:] /= rvec[:,None]
    grad[:,3:] /= rvec.size
    grad[:,:3] = -grad[:,3:]
    return sigma, grad.T

def bondangle(beads, indices, target):
    """Calculate the mean of the bond-angle over the replicas and
       return the value of the constraint function and its gradient.
    """
    grad = np.empty((beads.nbeads,len(indices)),float)
    for i,idx in enumerate(indices):
        grad[:,i] = dstrip(beads.q[:,idx])
    qXA = grad[:,3:6]-grad[:,:3]
    rXA = np.sqrt(np.sum(qXA**2, axis=-1))
    qXA /= rXA[:,None]
    qXB = grad[:,6:]-grad[:,:3]
    rXB = np.sqrt(np.sum(qXB**2, axis=-1))
    qXB /= rXB[:,None]
    ct = np.sum(qXA*qXB, axis=-1)
    st = np.sqrt(1-ct**2)
    sigma = np.arccos(ct).mean()
    sigma -= target
    grad[:,3:6] = (ct[:,None]*qXA-qXB)/(rXA*st*rXA.size)[:,None]
    grad[:,6:] = (ct[:,None]*qXB-qXA)/(rXB*st*rXB.size)[:,None]
    grad[:,:3] = -grad[:,3:6]
    grad[:,:3] -= grad[:,6:]

    return sigma, grad.T

def test_replicas():
    beads, replica_list = create_beads(local("test.ice_Ih.xyz"))
    assert beads.q[:,0] == pytest.approx(replica_list[0].q)
    assert beads.p[:,0] == pytest.approx(replica_list[0].p)
    assert beads.m[0] == pytest.approx(replica_list[0].m)
    beads.q[:,0] += 1
    assert not (beads.q[:,0] == pytest.approx(replica_list[0].q))
    replica_list[0].q += 1
    assert beads.q[:,0] == pytest.approx(replica_list[0].q)

def test_bond():
    beads, replicas = create_beads(local("test.ice_Ih.xyz"))
    dofs = list(range(6))
    bond_constraint = BondLength(dofs)
    bond_constraint.bind(replicas)
    sigma, grad = bondlength(beads, dofs, bond_constraint.targetval)
    assert bond_constraint.sigma == pytest.approx(sigma)
    assert bond_constraint.jac == pytest.approx(grad)
    # vary targetval
    new_target = 2.5
    bond_constraint.targetval = new_target
    sigma, grad = bondlength(beads, dofs, new_target)
    assert bond_constraint.sigma == pytest.approx(sigma)
    assert bond_constraint.jac == pytest.approx(grad)
    # vary configuration
    beads.q[0,dofs[2]] = 1.0
    sigma, grad = bondlength(beads, dofs, new_target)
    assert not (bond_constraint.sigma == pytest.approx(sigma))
    assert not (bond_constraint.jac == pytest.approx(grad))
    replicas[dofs[2]].q[0] = beads.q[0,dofs[2]]
    assert (bond_constraint.sigma == pytest.approx(sigma))
    assert (bond_constraint.jac == pytest.approx(grad))

def test_angle():
    beads, replicas = create_beads(local("test.ice_Ih.xyz"))
    dofs = list(range(9))
    bond_constraint = BondAngle(dofs)
    bond_constraint.bind(replicas)
    sigma, grad = bondangle(beads, dofs, bond_constraint.targetval)
    assert bond_constraint.sigma == pytest.approx(sigma)
    assert bond_constraint.jac == pytest.approx(grad)
    # vary targetval
    new_target = np.pi/4
    bond_constraint.targetval = new_target
    sigma, grad = bondangle(beads, dofs, new_target)
    assert bond_constraint.sigma == pytest.approx(sigma)
    assert bond_constraint.jac == pytest.approx(grad)
    # vary configuration
    beads.q[0,dofs[2]] = 1.0
    replicas[dofs[2]].q[0] = beads.q[0,dofs[2]]
    sigma, grad = bondangle(beads, dofs, new_target)
    assert bond_constraint.sigma == pytest.approx(sigma)
    assert bond_constraint.jac == pytest.approx(grad)