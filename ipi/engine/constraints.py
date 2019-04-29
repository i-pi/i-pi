"""Classes that deal with holonomic constraints.

Contains objects that return the values and the gradients of holonomic constraint
functions.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2019 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np

from ipi.utils.depend import dobject, dd, depend_array, depend_value, dstrip

__all__ = ['Replicas','HolonomicConstraint','BondLength','BondAngle']

class Replicas(dobject):

    """Storage of ring-polymer replica positions, masses and momenta
       for a single degree of freedom.

    Positions and momenta are stored as nbeads-sized contiguous arrays,
    and mass is stored as a scalar.

    Attributes:
       nbeads: The number of beads.

    Depend objects:
       p: An array giving the replica momenta for this DoF.
       q: An array giving the replica positions for this DoF.
       m: A scalar giving the replica masses
    """

    def __init__(self, nbeads, beads=None, idof=None):
        """Initialises Replicas.

        Args:
           nbeads: An integer giving the number of beads.
           beads: An instance of Beads from which to copy the initial values
                  of the positions, momenta and masses
           idof: An integer index of the degree of freedom from which to copy
                  these values.
        """
        self.nbeads = nbeads
        dself = dd(self)
        if beads is None:
            qtemp = np.zeros(nbeads, float)
            ptemp = np.zeros(nbeads, float)
            mtemp = 0.0
        else:
            if idof is None:
                raise ValueError("The index of the degree of freedom must be "+
                                 "specified when initialising Replicas from "+
                                 "beads.")
            qtemp = dstrip(beads.q[:,idof]).copy()
            ptemp = dstrip(beads.p[:,idof]).copy()
            mtemp = beads.m3[0,idof]

        dself.q = depend_array(name="q", value=qtemp)
        dself.p = depend_array(name="p", value=ptemp)
        dself.m = depend_value(name="m", value=mtemp)

    def __len__(self):
        """Length function.

        This is called whenever the standard function len(replicas) is used.

        Returns:
           The number of beads.
        """

        return self.nbeads

    def copy(self):
        """Creates a new Replicas object.

        Returns:
           A Replicas object with the same q, p, and m as the original.
        """

        newrep = Replicas(self.nbeads)
        newrep.q[:] = self.q
        newrep.p[:] = self.p
        newrep.m = self.m
        return newrep

class HolonomicConstraint(dobject):
    """Base holonomic constraints class.

    Specifies the standard methods and attributes common to all holonomic
    constraint classes.

    Attributes:
        nbeads: The number of beads in the set of constrained ring-polymer
                DoFs bound to the object
        ndof: The number of degrees of freedom involved in the constraint

    Depend objects:
        dofs: 1D array of integers indexing the degrees of freedom bound to
              the constraint
        targetval: The target value of the constraint function
        _q: List of Replicas' coordinates bound to the constraint
        _m: List of Replicas' masses bound to the constraint

    """
    # number of degrees of freedom expected by the constraint function
    _ndof = 0

    @classmethod
    def get_ndof(cls):
        """The number of degrees of freedom expected by the constraint function.
        """
        if (cls._ndof == 0):
            raise TypeError("Trying to get the number of degrees of freedom "+
                            "from base HolonomicConstraint class")
        return cls._ndof

    def __init__(self, dofs, val=None):
        """Initialise the holonomic constraint.
        """
        dself = dd(self)
        ndof = dself.get_ndof()
        dofs = np.asarray(dofs, dtype=int)
        if (dofs.ndim != 1):
            raise ValueError("Shape of constrained DoF group incompatible "+
                             "with "+dself.__class__.__name__)
        if (dofs.shape != (ndof,) and ndof != -1):
            raise ValueError("Size of constrained DoF group incompatible "+
                             "with "+dself.__class__.__name__)
        dself.dofs = depend_array(name="dofs",value=dofs)
        self.ndof = ndof
        if val is not None:
            if not np.isscalar(val):
                raise TypeError("Expecting a scalar for target constraint value")
        dself.targetval = depend_value(name="targetval")
        self.targetval = val

    def bind(self, replicas, **kwargs):
        """Bind the appropriate degrees of freedom to the holonomic constraint.
           Input:
               replicas(list): List of Replicas objects containing the ring-polymer
                               data grouped by degree of freedom.
        """

        self.nbeads = replicas[0].q.size
        dself = dd(self)
        dself._q = [replicas[i].q for i in self.dofs]
        dself._m = depend_array(name="_m",
                                value=np.array([replicas[i].m for i in self.dofs])
                                )

    def get_sigma(self):
        """Dummy constraint function calculator that does nothing."""
        pass

    def get_jac(self):
        """Dummy constraint Jacobian calculator that does nothing."""
        pass

#------------- Specialisations of the Holonomic Constraint class -----------#
class BondLength(HolonomicConstraint):
    """Constraint on the length of a bond between two atoms,
       averaged over the replicas.
    """

    # Six degrees of freedom --- three per atom
    _ndof = 6

    def bind(self, replicas, **kwargs):
        """Bind the appropriate degrees of freedom to the holonomic constraint.
        """

        super(BondLength, self).bind(replicas, **kwargs)
        dself = dd(self)
        dself._qAB = depend_array(name="_qAB", func=self.get_qAB,
                value=np.zeros((3,self.nbeads), float),
                dependencies=dself._q)
        dself._rAB = depend_array(name="_rAB", func=self.get_rAB,
                value=np.zeros(self.nbeads, float),
                dependencies=[dself._qAB])
        dself.sigma = depend_value(name="sigma", func=self.get_sigma,
                dependencies=dself._q+[dself.targetval])
        if self.targetval is None:
            # Set to the current value of the constraint function
            self.targetval = 0.0
            currentval = self.sigma
            self.targetval = currentval
        dself.jac = depend_value(name="jac", func=self.get_jac,
                                 dependencies=dself._q)

    def get_qAB(self):
        """Calculate the displacement vectors between the replicas of
           the two atoms.
        """
        return np.asarray([self._q[i+3]-self._q[i] for i in range(3)])

    def get_rAB(self):
        """Calculate the magnitude of the displacement vectors between the
           replicas of the two atoms.
        """
        return np.sqrt(np.sum(self._qAB**2, axis=0))

    def get_sigma(self):
        """Calculate the difference between the mean bond-length and the
           target value.
        """
        return np.mean(self._rAB) - self.targetval

    def get_jac(self):
        """Calculate the gradient of the constraint function.
        """
        jac = np.empty_like(dstrip(self._q))
        jac[3:] = self._qAB
        jac[3:] /= self._rAB
        jac[3:] /= self.nbeads
        jac[:3] = -jac[3:]
        return jac

class BondAngle(HolonomicConstraint):
    """Constraint on the A--X--B bond-angle, averaged over the replicas.
       The indices of the coordinates of the central atom X are supplied
       first, followes by A, then B.
    """

    # Nine degrees of freedom -- three per atom
    _ndof =  9

    def bind(self, replicas, **kwargs):
        """Bind the appropriate degrees of freedom to the holonomic constraint.
        """

        super(BondAngle,self).bind(replicas, **kwargs)
        dself = dd(self)
        dself._qXA = depend_array(name="_qXA", func=self.get_qXA,
                                  value=np.zeros((3,self.nbeads),float),
                                  dependencies=dself._q)
        dself._rXA = depend_array(name="_rXA", func=self.get_rXA,
                                  value=np.zeros(self.nbeads,float),
                                  dependencies=[dself._qXA])
        dself._unit_qXA = depend_array(name="_unit_qXA", func=self.get_unit_qXA,
                                  value=np.zeros((3,self.nbeads), float),
                                  dependencies=[dself._rXA])
        dself._qXB = depend_array(name="_qXB", func=self.get_qXB,
                                  value=np.zeros((3,self.nbeads), float),
                                  dependencies=dself._q)
        dself._rXB = depend_array(name="_rXB", func=self.get_rXB,
                                  value=np.zeros(self.nbeads, float),
                                  dependencies=[dself._qXB])
        dself._unit_qXB = depend_array(name="_unit_qXB", func=self.get_unit_qXB,
                                  value=np.zeros((3,self.nbeads), float),
                                  dependencies=[dself._rXB])
        dself._ct = depend_array(name="_ct", func=self.get_ct,
                                 value=np.zeros(self.nbeads, float),
                                 dependencies=[dself._unit_qXA, dself._unit_qXB])
        dself.sigma = depend_value(name="sigma", func=self.get_sigma,
                                   dependencies=dself._q+[dself.targetval])
        if self.targetval is None:
            # Set to initial bond angle
            self.targetval = 0.0
            currentval = self.sigma
            self.targetval = currentval
        dself.jac = depend_value(name="jac", func=self.get_jac,
                                 dependencies=dself._q)

    def get_qXA(self):
        """Calculate the displacement vectors between the replicas of
           the central atom and atom A
        """
        return np.asarray([self._q[i+3]-self._q[i] for i in range(3)])


    def get_rXA(self):
        """Calculate the magnitude of the displacement vectors between the
           replicas of the central atom and atom A
        """
        return np.sqrt(np.sum(self._qXA**2, axis=0))

    def get_unit_qXA(self):
        """Calculate the normalised displacement vectors between the replicas of
           the central atom and atom A
        """
        return self._qXA/self._rXA

    def get_qXB(self):
        """Calculate the displacement vectors between the replicas of
           the central atom and atom B
        """
        return np.asarray([self._q[i+6]-self._q[i] for i in range(3)])

    def get_rXB(self):
        """Calculate the magnitude of the displacement vectors between the
           replicas of the central atom and atom B
        """
        return np.sqrt(np.sum(self._qXB**2, axis=0))

    def get_unit_qXB(self):
        """Calculate the normalised displacement vectors between the replicas of
           the central atom and atom B
        """
        return self._qXB/self._rXB

    def get_ct(self):
        """Calculate the cosines of the bond-angles across the ring-polymer
           replicas.
        """
        return np.sum(self._unit_qXA*self._unit_qXB, axis=0)

    def get_sigma(self):
        """Calculate the difference between the mean angle and the
           target value.
        """
        return np.mean(np.arccos(self._ct)) - self.targetval

    def get_jac(self):
        """Calculate the gradient of the constraint function.
        """
        jac = np.zeros_like(self._q)
        st = np.sqrt(1-self._ct**2)
        jac[3:6] = (dstrip(self._ct)*dstrip(self._unit_qXA) -
                    dstrip(self._unit_qXB))
        jac[3:6] /= dstrip(self._rXA)
        jac[3:6] /= st
        jac[3:6] /= dstrip(self._rXA).size
        jac[6:] = (dstrip(self._ct)*dstrip(self._unit_qXB) -
                    dstrip(self._unit_qXA))
        jac[6:] /= dstrip(self._rXB)
        jac[6:] /= st
        jac[6:] /= dstrip(self._rXB).size
        jac[:3] -= jac[3:6]
        jac[:3] -= jac[6:]
        return jac

class EckartTrans(HolonomicConstraint):
    """Constraint on one of the components of the centre of mass.
    """

    # Number of DoFs determined upon initialisation
    _ndof = -1

    def __init__(self, dofs, coord, val=None):
        """Initialise the holonomic constraint.

           Args:
               dofs(list): integers indexing *all* the degrees of freedom of the
                           atoms subject to this Eckart constraint
               coord(str): 'x', 'y', 'z' -- specifies the component of the CoM
                           to be constrained
               val(float): the position at which the component is to be constrained.
        """
        q_str = coord.lower()
        if q_str=="x":
            idx = 0
        elif q_str=="y":
            idx = 1
        elif q_str=="z":
            idx = 2
        else:
            raise ValueError("Invalid coordinate specification supplied to "+
                             self.__class__.__name__)
        super(EckartTrans,self).__init__(dofs[idx::3], val)

    def bind(self, replicas, **kwargs):
        """Bind the appropriate coordinates to the constraint.
        """

        super(EckartTrans, self).bind(replicas, **kwargs)
        dself = dd(self)
        dself._qc = depend_array(name="_qc", func=self.get_qc,
                                 value=np.zeros(self.ndof, float),
                                 dependencies=dself._q)
        dself._mrel = depend_array(name="_mrel", func=self.get_mrel,
                                 value=np.zeros(self.ndof, float),
                                 dependencies=dself._m)
        dself.sigma = depend_value(name="sigma", func=self.get_sigma,
                                   dependencies=[dself._qc, dself._mrel,
                                                 dself.targetval])
        if self.targetval is None:
            # Set to the current position of the CoM
            self.targetval = 0.0
            currentval = self.sigma
            self.targetval = currentval
        dself.jac = depend_value(name="jac", func=self.get_jac,
                                 dependencies=[dself._mrel])

    def get_qc(self):
        """Calculate the centroids of the bead DoFs.
        """

        return np.asarray([np.mean(q) for q in self._q])

    def get_mrel(self):
        """Return the masses of the centroids, divided by the total mass.
        """

        ans = np.asarray([m[0] for m in self._m])
        ans /= ans.sum()
        return ans

    def get_sigma(self):
        return np.dot(self._mrel, self._qc) - self.targetval

    def get_jac(self):
        jac = np.empty_like(dstrip(self._q))
        jac[...] = self._mrel[:,None]
        jac /= self.nbeads
        return jac

class EckartTransX(EckartTrans):
    """Constraint on the x-component of the CoM of a group of atoms.
    """
    def __init__(self, dofs, val=None):
        """Initialise the holonomic constraint.

           Args:
               dofs(list): integers indexing *all* the degrees of freedom of the
                           atoms subject to this Eckart constraint
               val(float): the position at which the component is to be constrained.
        """
        super(EckartTransX,self).__init__(dofs, "x", val)

class EckartTransY(EckartTrans):
    """Constraint on the y-component of the CoM of a group of atoms.
    """
    def __init__(self, dofs, val=None):
        """Initialise the holonomic constraint.

           Args:
               dofs(list): integers indexing *all* the degrees of freedom of the
                           atoms subject to this Eckart constraint
               val(float): the position at which the component is to be constrained.
        """
        super(EckartTransY,self).__init__(dofs, "y", val)

class EckartTransZ(EckartTrans):
    """Constraint on the z-component of the CoM of a group of atoms.
    """
    def __init__(self, dofs, val=None):
        """Initialise the holonomic constraint.

           Args:
               dofs(list): integers indexing *all* the degrees of freedom of the
                           atoms subject to this Eckart constraint
               val(float): the position at which the component is to be constrained.
        """
        super(EckartTransZ,self).__init__(dofs, "z", val)

class EckartRot(HolonomicConstraint):
    """One of the components of the Eckart rotational constraint.

       NOTE: in this definition the usual sum over cross-products is divided by
             the total mass of the system.
    """

    # Number of DoFs determined upon initialisation
    _ndof = -1

    def __init__(self, dofs, coord, val=None):
        """Initialise the holonomic constraint.

           Args:
               dofs(list): integers indexing *all* the degrees of freedom of the
                           atoms subject to this Eckart constraint
               coord(str): 'x', 'y', 'z' -- specifies the component of the
                           angular-momentum-like quantity to be constrained
               val(array-like): this argument is not used
        """
        q_str = coord.lower()
        if q_str=="x":
            idces = (1,2)
        elif q_str=="y":
            idces = (2,0)
        elif q_str=="z":
            idces = (0,1)
        else:
            raise ValueError("Invalid coordinate specification supplied to "+
                             self.__class__.__name__)
        super(EckartRot,self).__init__(dofs[idces[0]::3]+dofs[idces[1]::3], 0.0)

    def bind(self, replicas, **kwargs):
        """Bind the appropriate coordinates to the constraints.
          Args:
              replicas(list): List of Replicas
          **kwargs:
              ref(array-like): Reference configuration for the constraint; if
                               absent, taken to be the centroid configuration

        """
        super(EckartRot, self).bind(replicas, **kwargs)
        dself = dd(self)
        dself._mrel = depend_array(name="_mrel", func=self.get_mrel,
                                 value=np.zeros((2,self.ndof//2), float),
                                 dependencies=dself._m)
        if "ref" in kwargs:
            lref = np.asarray(kwargs["ref"]).flatten()
            ref = np.asarray(
                    [ lref[i] for i in self.dofs ]).reshape((2, self.ndof//2))
        else:
            # Use the centroids
            ref = np.asarray(
                    [ np.mean(q) for q in self._q] ).reshape((2, self.ndof//2))
        dself._ref = depend_array(name="_ref", value=ref)
        #NOTE: this is the mass-weighted configuration *in its CoM*
        dself._mref = depend_array(name="_mref", func=self.get_mref,
                                   value=dstrip(self._mrel)*dstrip(self._ref),
                                   dependencies=[dself._mrel, dself._ref])
        # Displacements between centroids and reference
        dself._delqc = depend_array(name="_delqc", func=self.get_delqc,
                                 value=np.zeros((2,self.ndof//2), float),
                                 dependencies=dself._q+[dself._ref])
        dself.sigma = depend_value(name="sigma", func=self.get_sigma,
                                   dependencies=[dself._delqc, dself._mref])
        dself.jac = depend_value(name="jac", func=self.get_jac,
                                 dependencies=[dself._mref])

    def get_mrel(self):
        """Return the masses of the centroids, divided by the total mass.
        """

        ans = np.asarray([m[0] for m in self._m]).reshape((2,self.ndof//2))
        ans /= ans.sum(axis=-1)[:,None]
        return ans

    def get_mref(self):
        """Return the coordinates of the reference configuration in its CoM frame,
           weighted by self._mrel
        """
        ans = self._mrel*self._ref
        CoM = dstrip(ans).sum(axis=-1)
        ans[...] = dstrip(self._ref)-CoM[:,None]
        ans *= dstrip(self._mrel)
        return ans

    def get_delqc(self):
        """Return the displacement between the centroid coordinates and the
           reference configuration (not mass-weighted, lab-frame).
        """
        delqc = np.asarray(
                    [ np.mean(q) for q in self._q] ).reshape((2, self.ndof//2))
        delqc -= self._ref
        return delqc

    def get_sigma(self):
        return np.sum(self._mref[0]*self._delqc[1] - self._mref[1]*self._delqc[0])

    def get_jac(self):
        jac = np.empty_like(dstrip(self._q))
        jac[:self.ndof//2,:] = -self._mref[:,1,None]
        jac[self.ndof//2:,:] = self._mref[:,0,None]
        jac /= self.nbeads
        return jac

class EckartRotX(EckartRot):
    """Constraint on the x-component of the Eckart "angular momentum"
    """
    def __init__(self, dofs, val=None):
        """Initialise the holonomic constraint.

           Args:
               dofs(list): integers indexing *all* the degrees of freedom of the
                           atoms subject to this Eckart constraint
               val(float): the position at which the component is to be constrained.
        """
        super(EckartRotX,self).__init__(dofs, "x", val)

class EckartRotY(EckartRot):
    """Constraint on the y-component of the Eckart "angular momentum"
    """
    def __init__(self, dofs, val=None):
        """Initialise the holonomic constraint.

           Args:
               dofs(list): integers indexing *all* the degrees of freedom of the
                           atoms subject to this Eckart constraint
               val(float): the position at which the component is to be constrained.
        """
        super(EckartRotY,self).__init__(dofs, "y", val)

class EckartRotZ(EckartRot):
    """Constraint on the z-component of the Eckart "angular momentum".
    """
    def __init__(self, dofs, val=None):
        """Initialise the holonomic constraint.

           Args:
               dofs(list): integers indexing *all* the degrees of freedom of the
                           atoms subject to this Eckart constraint
               val(float): the position at which the component is to be constrained.
        """
        super(EckartRotZ,self).__init__(dofs, "z", val)
