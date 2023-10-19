"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np

from ipi.utils.messages import verbosity,warning
from ipi.utils.depend import *
from ipi.utils.units import Constants

__all__ = ['DummyIntegrator', \
           'NVEIntegrator', 'NVTIntegrator', 'NPTIntegrator', \
           'EDANVEIntegrator' , 'EDANVTIntegrator', \
           'SCIntegrator', 'SCNPTIntegrator', 'NSTIntegrator', 'NVTCCIntegrator']

class DummyIntegrator(dobject):
    """No-op integrator for (PI)MD"""

    def __init__(self):
        #super(TimeDependentIntegrator,self).__init__()
        pass

    def get_qdt(self):
        return self.dt * 0.5 / self.inmts

    def get_pdt(self):
        dtl = 1.0 / self.nmts
        for i in range(1, len(dtl)):
            dtl[i] *= dtl[i - 1]
        dtl *= self.dt * 0.5
        return dtl

    def get_tdt(self):
        if self.splitting == "obabo":
            return self.dt * 0.5
        elif self.splitting == "baoab":
            return self.dt
        else:
            raise ValueError(
                "Invalid splitting requested. Only OBABO and BAOAB are supported."
            )

    def bind(self, motion):
        """Reference all the variables for simpler access."""
        
        self.beads      = motion.beads
        self.bias       = motion.ensemble.bias
        self.ensemble   = motion.ensemble
        self.forces     = motion.forces
        self.prng       = motion.prng
        self.nm         = motion.nm
        self.thermostat = motion.thermostat
        self.barostat   = motion.barostat
        self.fixcom     = motion.fixcom
        self.fixatoms   = motion.fixatoms
        self.enstype    = motion.enstype

        dself = dd(self)
        dmotion = dd(motion)

        # no need to dpipe these are really just references
        dself.splitting = dmotion.splitting
        dself.dt = dmotion.dt
        dself.nmts = dmotion.nmts

        # total number of iteration in the inner-most MTS loop
        dself.inmts = depend_value(name="inmts", func=lambda: np.prod(self.nmts))
        dself.nmtslevels = depend_value(name="nmtslevels", func=lambda: len(self.nmts))
        # these are the time steps to be used for the different parts of the integrator
        dself.qdt = depend_value(
            name="qdt",
            func=self.get_qdt,
            dependencies=[dself.splitting, dself.dt, dself.inmts],
        )  # positions
        dself.pdt = depend_array(
            name="pdt",
            func=self.get_pdt,
            value=np.zeros(len(self.nmts)),
            dependencies=[dself.splitting, dself.dt, dself.nmts],
        )  # momenta
        dself.tdt = depend_value(
            name="tdt",
            func=self.get_tdt,
            dependencies=[dself.splitting, dself.dt, dself.nmts],
        )  # thermostat

        dpipe(dself.qdt, dd(self.nm).dt)
        dpipe(dself.dt , dd(self.barostat).dt)
        dpipe(dself.qdt, dd(self.barostat).qdt)
        dpipe(dself.pdt, dd(self.barostat).pdt)
        dpipe(dself.tdt, dd(self.barostat).tdt)
        dpipe(dself.tdt, dd(self.thermostat).dt)

        if motion.enstype == "sc" or motion.enstype == "scnpt":
            # coefficients to get the (baseline) trotter to sc conversion
            self.coeffsc = np.ones((self.beads.nbeads, 3 * self.beads.natoms), float)
            self.coeffsc[::2] /= -3.0
            self.coeffsc[1::2] /= 3.0

    def pstep(self):
        """Dummy momenta propagator which does nothing."""
        pass

    def qcstep(self):
        """Dummy centroid position propagator which does nothing."""
        pass

    def step(self, step=None):
        """Dummy simulation time step which does nothing."""
        pass

    def pconstraints(self):
        """This removes the centre of mass contribution to the kinetic energy.

        Calculates the centre of mass momenta, then removes the mass weighted
        contribution from each atom. If the ensemble defines a thermostat, then
        the contribution to the conserved quantity due to this subtraction is
        added to the thermostat heat energy, as it is assumed that the centre of
        mass motion is due to the thermostat.

        If there is a choice of thermostats, the thermostat
        connected to the centroid is chosen.
        """

        if self.fixcom:
            na3 = self.beads.natoms * 3
            nb = self.beads.nbeads
            p = dstrip(self.beads.p)
            m = dstrip(self.beads.m3)[:, 0:na3:3]
            M = self.beads[0].M
            Mnb = M * nb

            dens = 0
            for i in range(3):
                pcom = p[:, i:na3:3].sum()
                dens += pcom**2
                pcom /= Mnb
                self.beads.p[:, i:na3:3] -= m * pcom

            self.ensemble.eens += dens * 0.5 / Mnb

        if len(self.fixatoms) > 0:
            for bp in self.beads.p:
                m = dstrip(self.beads.m)
                self.ensemble.eens += 0.5 * np.dot(
                    bp[self.fixatoms * 3], bp[self.fixatoms * 3] / m[self.fixatoms]
                )
                self.ensemble.eens += 0.5 * np.dot(
                    bp[self.fixatoms * 3 + 1],
                    bp[self.fixatoms * 3 + 1] / m[self.fixatoms],
                )
                self.ensemble.eens += 0.5 * np.dot(
                    bp[self.fixatoms * 3 + 2],
                    bp[self.fixatoms * 3 + 2] / m[self.fixatoms],
                )
                bp[self.fixatoms * 3] = 0.0
                bp[self.fixatoms * 3 + 1] = 0.0
                bp[self.fixatoms * 3 + 2] = 0.0

class NVEIntegrator(DummyIntegrator):

    """Integrator object for constant energy simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant energy ensemble. Note that a temperature of some kind must be
    defined so that the spring potential can be calculated.

    Attributes:
        ptime: The time taken in updating the velocities.
        qtime: The time taken in updating the positions.
        ttime: The time taken in applying the thermostat steps.

    Depend objects:
        econs: Conserved energy quantity. Depends on the bead kinetic and
            potential energy, and the spring potential energy.
    """

    def pstep(self, level=0):
        """Velocity Verlet momentum propagator."""

        # print(" # NVEIntegrator.pstep")
        # halfdt/alpha
        a = self.forces.forces_mts(level) * self.pdt[level] 
        self.beads.p += a
        if level == 0:  # adds bias in the outer loop
            self.beads.p += dstrip(self.bias.f) * self.pdt[level]

    def qcstep(self):
        """Velocity Verlet centroid position propagator."""
        # dt/inmts
        self.nm.qnm[0, :] += (dstrip(self.nm.pnm)[0, :] / dstrip(self.beads.m3)[0] * self.qdt)

    # now the idea is that for BAOAB the MTS should work as follows:
    # take the BAB MTS, and insert the O in the very middle. This might imply breaking a A step in two, e.g. one could have
    # Bbabb(a/2) O (a/2)bbabB
    def mtsprop_ba(self, index):
        """Recursive MTS step"""

        mk = int(self.nmts[index] / 2)

        for i in range(mk):  # do nmts/2 full sub-steps

            self.pstep(index)
            self.pconstraints()
            if index == self.nmtslevels - 1:
                # call Q propagation for dt/alpha at the inner step
                self.qcstep()
                self.nm.free_qstep()
                self.qcstep()
                self.nm.free_qstep()

            else:
                self.mtsprop(index + 1)

            self.pstep(index)
            self.pconstraints()

        if self.nmts[index] % 2 == 1:
            # propagate p for dt/2alpha with force at level index
            self.pstep(index)
            self.pconstraints()
            if index == self.nmtslevels - 1:
                # call Q propagation for dt/alpha at the inner step
                self.qcstep()
                self.nm.free_qstep()
            else:
                self.mtsprop_ba(index + 1)

    def mtsprop_ab(self, index):
        """Recursive MTS step"""

        if self.nmts[index] % 2 == 1:
            if index == self.nmtslevels - 1:
                # call Q propagation for dt/alpha at the inner step
                self.qcstep()
                self.nm.free_qstep()
            else:
                self.mtsprop_ab(index + 1)

            # propagate p for dt/2alpha with force at level index
            self.pstep(index)
            self.pconstraints()

        for i in range(int(self.nmts[index] / 2)):  # do nmts/2 full sub-steps
            self.pstep(index)
            self.pconstraints()
            if index == self.nmtslevels - 1:
                # call Q propagation for dt/alpha at the inner step
                self.qcstep()
                self.nm.free_qstep()
                self.qcstep()
                self.nm.free_qstep()
            else:
                self.mtsprop(index + 1)

            self.pstep(index)
            self.pconstraints()

    def mtsprop(self, index):
        # just calls the two pieces together
        self.mtsprop_ba(index) # pstep, qcstep
        self.mtsprop_ab(index) # qcstep, pstep

    def step(self, step=None):
        """Does one simulation time step."""
        # print(" # NVEIntegrator.step")
        self.mtsprop(0)

class NVTIntegrator(NVEIntegrator):

    """Integrator object for constant temperature simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant temperature ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant.

    Attributes:
        thermostat: A thermostat object to keep the temperature constant.
    """

    def tstep(self):
        """Velocity Verlet thermostat step"""

        self.thermostat.step()

    def step(self, step=None):
        """Does one simulation time step."""

        if self.splitting == "obabo":
            # thermostat is applied for dt/2
            self.tstep()
            self.pconstraints()

            # forces are integerated for dt with MTS.
            self.mtsprop(0)

            # thermostat is applied for dt/2
            self.tstep()
            self.pconstraints()

        elif self.splitting == "baoab":

            self.mtsprop_ba(0)
            # thermostat is applied for dt
            self.tstep()
            self.pconstraints()
            self.mtsprop_ab(0)

class NPTIntegrator(NVTIntegrator):

    """Integrator object for constant pressure simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant pressure ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant, and a barostat to keep the
    pressure constant.
    """

    # should be enough to redefine these functions, and the step() from NVTIntegrator should do the trick

    def pstep(self, level=0):
        """Velocity Verlet monemtum propagator."""

        if np.array_equiv(self.forces.vir, np.zeros(len(self.forces.vir))):
            raise ValueError(
                "Seems like no stress tensor was computed by the client. Stopping barostat!"
            )
        self.barostat.pstep(level)
        super(NPTIntegrator, self).pstep(level)
        # self.pconstraints()

    def qcstep(self):
        """Velocity Verlet centroid position propagator."""

        self.barostat.qcstep()

    def tstep(self):
        """Velocity Verlet thermostat step"""

        self.thermostat.step()
        self.barostat.thermostat.step()
        # self.pconstraints()

class EDAIntegrator(DummyIntegrator):
    """Integrator object for simulations using the electric dipole approximation when an external electric field is applied.

    The electric field has to be assumed:
        1) homogeneous all over the system
        2) slowly varying.
    The first assumption need to be done in order for the Electric Dipole Approximation to be valid.
    The second one in stead is neeed in order to DFT calculations to be reliable, 
    since they are performed without any external electric field applied 
    (to keep the KS equations simple(r) and computationally affordable).

    """

    # author: Elia Stocco
    # motivation: deal with time-dependent external potential, in particular the EDA potential (have a look at EDAIntegrator)
    # information: not included in __all__

    # the conversion factor for the polarization from C/m^2 to atomic unit
    # 1e    = 1.602176634×10−19 C
    # 1Bohr = 5.291772109×10−11 m
    # C/m^2 = (5.291772109×10−11)^2/1.602176634×10−19 = 1.7478E-2 
    # from_Cm2_to_atomic = 1.7478E-2 
    # from_atomic_to_Cm2 = 57.214766

    # threshold for the comparison between the ionic polarization computed by the driven and the one computed by i-pi
    # pay attention that they could be numerically different but anyway equivalent due to a different branch mapping
    # thr_pol_comparison  = 1.0

    # if true, _check = _check_dipole, otherwise a check of the polarizations values is performed too
    # _check_flag = False

    def __init__(self):
        super(EDAIntegrator,self).__init__()

    def bind(self,motion):
        """bind variables"""
        super(EDAIntegrator,self).bind(motion) 

        dself = dd(self)

        dep = [dd(self.ensemble).time,dd(self.ensemble.eda).cptime,dd(self.ensemble.eda).bec,dd(self.ensemble.eda).Efield]
        dself.EDAforces  = depend_array(name="forces"  , func=self._forces              ,\
                                        value=np.zeros((self.beads.nbeads,self.beads.natoms*3)),dependencies=dep)
        dself.EDAxforces = depend_array(name="xforces" , func=lambda:self._forces_component("x") ,\
                                        value=np.zeros((self.beads.nbeads,self.beads.natoms))  ,dependencies=[dd(self).EDAforces])
        dself.EDAyforces = depend_array(name="yforces" , func=lambda:self._forces_component("y") ,\
                                        value=np.zeros((self.beads.nbeads,self.beads.natoms))  ,dependencies=[dd(self).EDAforces])
        dself.EDAzforces = depend_array(name="zforces" , func=lambda:self._forces_component("z") ,\
                                        value=np.zeros((self.beads.nbeads,self.beads.natoms))  ,dependencies=[dd(self).EDAforces])
        
        dself.dt = depend_value(name="dt", value=0.0)
        dpipe(dfrom=dd(motion).dt,dto=dself.dt)

        pass

    def pstep(self, level=0):
        """Velocity Verlet momentum propagator."""
        self.beads.p += self.EDAforces * self.pdt[level]
        pass

    def _forces(self):
        """Compute the EDA contribution to the forces due to the polarization"""
        # ES: if 'nbeads' > 1, then we will have to fix something here
        Z = self.ensemble.eda.bec
        E = dd(self.ensemble.eda).Efield(self.ensemble.eda.time)
        forces = Constants.e * Z @ E
        return forces 
        # return forces.reshape((self.beads.nbeads,-1))
    
    def _forces_component(self,xyz):
        """Get a component of the EDA contribution to the forces"""
        if xyz not in ["x","y","z"]:
            raise ValueError("'xyz' can assume only the values ['x','y','z']")
        if self.beads.nbeads != 1 :
            raise ValueError("'_forces_component' implemented only for nbeads=1")
        
        f = self.EDAforces[0].reshape((-1,3))
        if xyz == "x":
            return f[:,0]
        elif xyz == "y":
            return f[:,1]
        elif xyz == "z":
            return f[:,2]
        else :
            raise ValueError("'xyz' can assume only the values ['x','y','z']")

class EDANVEIntegrator(EDAIntegrator,NVEIntegrator):
    """Integrator object for simulations with constant number of particles, energy, and volume, using the electric dipole approximation when an external electric field is applied.

    The electric field has to be assumed:
        1) homogeneous all over the system
        2) slowly varying.
    The first assumption need to be done in order for the Electric Dipole Approximation to be valid.
    The second one in stead is neeed in order to DFT calculations to be reliable, 
    since they are performed without any external electric field applied 
    (to keep the KS equations simple(r) and computationally affordable).

    """

    # author: Elia Stocco
    # motivation: deal with time-dependent external potential, in particular the EDA potential (have a look at EDAIntegrator)

    # order = True

    # def pstep(self,level):
    #     if self.order :
    #         EDAIntegrator.pstep(self,level)
    #         NVEIntegrator.pstep(self,level)
    #         self.order = False
    #     else :
    #         NVEIntegrator.pstep(self,level) # the drive is called here
    #         EDAIntegrator.pstep(self,level)
    #         self.order = True
    #     #self.cptime += self.pdt[level] # update cptime
    #     pass

    def pstep(self,level):
        NVEIntegrator.pstep(self,level) # the driver is called here
        EDAIntegrator.pstep(self,level)
        # self.ensemble.eda.tacc   += dd(self.ensemble.eda).TderEDAenergy(self.ensemble.eda.cptime) * self.pdt[level] # wrong
        # self.ensemble.eda.cptime += self.pdt[level] # update cptime
        pass

    def step(self,step=None):
        super(EDANVEIntegrator,self).step(step)
        self.ensemble.eda.tacc += dd(self.ensemble.eda).TderEDAenergy(self.ensemble.eda.time) * self.dt

class EDANVTIntegrator(EDAIntegrator,NVTIntegrator):
    """Integrator object for simulations with constant number of particles, temperature, and volume, using the electric dipole approximation when an external electric field is applied.

    The electric field has to be assumed:
        1) homogeneous all over the system
        2) slowly varying.
    The first assumption need to be done in order for the Electric Dipole Approximation to be valid.
    The second one in stead is neeed in order to DFT calculations to be reliable, 
    since they are performed without any external electric field applied 
    (to keep the KS equations simple(r) and computationally affordable).

    """

    # author: Elia Stocco
    # motivation: deal with time-dependent external potential, in particular the EDA potential (have a look at EDAIntegrator)

    order = True

    def pstep(self,level):
        if self.order :
            EDAIntegrator.pstep(self,level)
            NVTIntegrator.pstep(self,level)
            self.order = False
        else :
            NVTIntegrator.pstep(self,level) # the drive is called here
            EDAIntegrator.pstep(self,level)
            self.order = True
        #self.cptime += self.pdt[level] # update cptime
        pass

class NVTCCIntegrator(NVTIntegrator):
    """Integrator object for constant temperature simulations with constrained centroid.

    Has the relevant conserved quantity and normal mode propagator for the
    constant temperature ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant.

    Attributes:
        thermostat: A thermostat object to keep the temperature constant.
    """

    def pstep(self):
        """Velocity Verlet momenta propagator."""

        # propagates in NM coordinates
        self.nm.pnm += dstrip(self.nm.fnm) * (self.dt * 0.5)
        # self.beads.p += dstrip(self.forces.f)*(self.dt*0.5)
        # also adds the bias force
        # self.beads.p += dstrip(self.bias.f)*(self.dt*0.5)

    def step(self, step=None):
        """Does one simulation time step."""

        self.thermostat.step()
        self.pconstraints()
        # NB we only have to take into account the energy balance of zeroing centroid velocity when we had added energy through the thermostat
        self.ensemble.eens += 0.5 * np.dot(
            self.nm.pnm[0], self.nm.pnm[0] / self.nm.dynm3[0]
        )
        self.nm.pnm[0, :] = 0.0

        self.pstep()
        self.nm.pnm[0, :] = 0.0
        self.pconstraints()

        # self.qcstep() # for the moment I just avoid doing the centroid step.
        self.nm.free_qstep()

        self.pstep()
        self.nm.pnm[0, :] = 0.0
        self.pconstraints()

        self.thermostat.step()
        self.ensemble.eens += 0.5 * np.dot(
            self.nm.pnm[0], self.nm.pnm[0] / self.nm.dynm3[0]
        )
        self.nm.pnm[0, :] = 0.0
        self.pconstraints()

class NSTIntegrator(NPTIntegrator):

    """Ensemble object for constant pressure simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant pressure ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant, and a barostat to keep the
    pressure constant.

    Attributes:
    barostat: A barostat object to keep the pressure constant.

    Depend objects:
    econs: Conserved energy quantity. Depends on the bead and cell kinetic
    and potential energy, the spring potential energy, the heat
    transferred to the beads and cell thermostat, the temperature and
    the cell volume.
    pext: External pressure.
    """

class SCIntegrator(NVTIntegrator):
    """Integrator object for constant temperature simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant temperature ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant.

    Attributes:
        thermostat: A thermostat object to keep the temperature constant.

    Depend objects:
        econs: Conserved energy quantity. Depends on the bead kinetic and
            potential energy, the spring potential energy and the heat
            transferred to the thermostat.
    """

    def bind(self, mover):
        """Binds ensemble beads, cell, bforce, bbias and prng to the dynamics.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the ensemble.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
        beads: The beads object from whcih the bead positions are taken.
        nm: A normal modes object used to do the normal modes transformation.
        cell: The cell object from which the system box is taken.
        bforce: The forcefield object from which the force and virial are
            taken.
        prng: The random number generator object which controls random number
            generation.
        """

        super(SCIntegrator, self).bind(mover)
        self.ensemble.add_econs(dd(self.forces).potsc)
        self.ensemble.add_xlpot(dd(self.forces).potsc)

    def pstep(self, level=0):
        """Velocity Verlet monemtum propagator."""

        if level == 0:
            # bias goes in the outer loop
            self.beads.p += dstrip(self.bias.f) * self.pdt[level]
        # just integrate the Trotter force scaled with the SC coefficients, which is a cheap approx to the SC force
        self.beads.p += (
            self.forces.forces_mts(level)
            * (1.0 + self.forces.coeffsc_part_1)
            * self.pdt[level]
        )

    def step(self, step=None):

        # the |f|^2 term is considered to be slowest (for large enough P) and is integrated outside everything.
        # if nmts is not specified, this is just the same as doing the full SC integration

        if self.splitting == "obabo":
            # thermostat is applied for dt/2
            self.tstep()
            self.pconstraints()

            # forces are integerated for dt with MTS.
            self.beads.p += dstrip(self.forces.fsc_part_2) * self.dt * 0.5
            self.mtsprop(0)
            self.beads.p += dstrip(self.forces.fsc_part_2) * self.dt * 0.5

            # thermostat is applied for dt/2
            self.tstep()
            self.pconstraints()

        elif self.splitting == "baoab":

            self.beads.p += dstrip(self.forces.fsc_part_2) * self.dt * 0.5
            self.mtsprop_ba(0)
            # thermostat is applied for dt
            self.tstep()
            self.pconstraints()
            self.mtsprop_ab(0)
            self.beads.p += dstrip(self.forces.fsc_part_2) * self.dt * 0.5

class SCNPTIntegrator(SCIntegrator):
    """Integrator object for constant pressure Suzuki-Chin simulations.

    Has the relevant conserved quantity and normal mode propagator for the
    constant pressure ensemble. Contains a thermostat object containing the
    algorithms to keep the temperature constant, and a barostat to keep the
    pressure constant.
    """

    # should be enough to redefine these functions, and the step() from NVTIntegrator should do the trick
    def pstep(self, level=0):
        """Velocity Verlet monemtum propagator."""

        if np.array_equiv(self.forces.vir, np.zeros(len(self.forces.vir))):
            raise ValueError(
                "Seems like no stress tensor was computed by the client. Stopping barostat!"
            )

        self.barostat.pstep(level)
        super(SCNPTIntegrator, self).pstep(level)

    def qcstep(self):
        """Velocity Verlet centroid position propagator."""

        self.barostat.qcstep()

    def tstep(self):
        """Velocity Verlet thermostat step"""

        self.thermostat.step()
        self.barostat.thermostat.step()

    def step(self, step=None):

        # the |f|^2 term is considered to be slowest (for large enough P) and is integrated outside everything.
        # if nmts is not specified, this is just the same as doing the full SC integration

        if self.splitting == "obabo":
            # thermostat is applied for dt/2
            self.tstep()
            self.pconstraints()

            # forces are integerated for dt with MTS.
            self.barostat.pscstep()
            self.beads.p += dstrip(self.forces.fsc_part_2) * self.dt * 0.5
            self.mtsprop(0)
            self.barostat.pscstep()
            self.beads.p += dstrip(self.forces.fsc_part_2) * self.dt * 0.5

            # thermostat is applied for dt/2
            self.tstep()
            self.pconstraints()

        elif self.splitting == "baoab":

            self.barostat.pscstep()
            self.beads.p += dstrip(self.forces.fsc_part_2) * self.dt * 0.5
            self.mtsprop_ba(0)
            # thermostat is applied for dt
            self.tstep()
            self.pconstraints()
            self.mtsprop_ab(0)
            self.barostat.pscstep()
            self.beads.p += dstrip(self.forces.fsc_part_2) * self.dt * 0.5
