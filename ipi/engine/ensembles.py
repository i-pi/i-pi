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

from ipi.utils.messages import warning,verbosity
from ipi.utils.depend import dd
from ipi.utils.units import Constants, UnitMap
from ipi.engine.thermostats import *
from ipi.engine.barostats import *
from ipi.engine.motion.alchemy import *
from ipi.engine.forces import Forces, ScaledForceComponent
import re

__all__ = ["Ensemble", "ensemble_swap"]

# IMPORTANT - THIS MUST BE KEPT UP-TO-DATE WHEN THE ENSEMBLE CLASS IS CHANGED


def ensemble_swap(ens1, ens2):
    """Swaps the definitions of the two ensembles, by
    exchanging all of the inner properties."""

    if ens1.temp != ens2.temp:
        ens1.temp, ens2.temp = ens2.temp, ens1.temp
    if ens1.pext != ens2.pext:
        ens1.pext, ens2.pext = ens2.pext, ens1.pext
    if np.linalg.norm(ens1.stressext - ens2.stressext) > 1e-10:
        tmp = dstrip(ens1.stressext).copy()
        ens1.stressext[:] = ens2.stressext
        ens2.stressext[:] = tmp
    if len(ens1.bweights) != len(ens2.bweights):
        raise ValueError(
            "Cannot exchange ensembles that have different numbers of bias components"
        )
    if len(ens1.hweights) != len(ens2.hweights):
        raise ValueError(
            "Cannot exchange ensembles that are described by different forces"
        )
    if not np.array_equal(ens1.bweights, ens2.bweights):
        ens1.bweights, ens2.bweights = (
            dstrip(ens2.bweights).copy(),
            dstrip(ens1.bweights).copy(),
        )
    if not np.array_equal(ens1.hweights, ens2.hweights):
        ens1.hweights, ens2.hweights = (
            dstrip(ens2.hweights).copy(),
            dstrip(ens1.hweights).copy(),
        )


class Ensemble(dobject):

    """Base ensemble class.

    Defines the thermodynamic state of the system.

    Depend objects:
        temp: The system's temperature.
        pext: The systems's pressure
        stressext: The system's stress tensor
        bias: Explicit bias forces
    """

    def __init__(
        self,
        eens=0.0,
        econs=0.0,
        temp=None,
        pext=None,
        stressext=None,
        bcomponents=None,
        bweights=None,
        hweights=None,
        time=0.0,
        Eamp=None,
        Efreq=None,
        Ephase=None,
        Epeak=None,
        Esigma=None,
        bec=None,
        cdip=True,
        # tacc=0.0,
        cbec=False,
    ):
        """Initialises Ensemble.

        Args:
            temp: The temperature.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """
        dself = dd(self)

        dself.temp = depend_value(name="temp")
        if temp is not None:
            self.temp = temp
        else:
            self.temp = -1.0

        dself.stressext = depend_array(name="stressext", value=np.zeros((3, 3), float))
        if stressext is not None:
            self.stressext = np.reshape(np.asarray(stressext), (3, 3))
        else:
            self.stressext = -1.0

        dself.pext = depend_value(name="pext")
        if pext is not None:
            self.pext = pext
        else:
            self.pext = -1.0

        dself.eens = depend_value(name="eens")
        if eens is not None:
            self.eens = eens
        else:
            self.eens = 0.0

        # the bias force contains two bits: explicit biases (that are meant to represent non-physical external biasing potentials)
        # and hamiltonian weights (that will act by scaling different physical components of the force). Both are bound as components
        # of the "bias force" evaluator, but their meaning (and the wiring further down in bind()) differ.

        # these are the additional bias components
        if bcomponents is None:
            bcomponents = []
        self.bcomp = bcomponents
        self.bias = Forces()

        # and their weights
        if bweights is None or len(bweights) == 0:
            bweights = np.ones(len(self.bcomp))

        dself.bweights = depend_array(name="bweights", value=np.asarray(bweights))

        # weights of the Hamiltonian scaling
        if hweights is None:
            hweights = np.ones(0)
        self.hweights = np.asarray(hweights)

        # ES

        if Epeak is not None and Epeak < 0 :
            raise ValueError("Epeak < 0: the peak of the external electric field can only be positive")    
        if Esigma is not None and Esigma < 0 :
            raise ValueError("Esigma < 0: the standard deviation of the gaussian envelope function of the external electric field has to be positive") 

        # Internal time counter
        dself.time = depend_value(name="time",value=time)

        #
        dself.Eamp   = depend_array(name="Eamp"  ,value=Eamp   if Eamp   is not None else np.zeros(3))
        dself.Efreq  = depend_value(name="Efreq" ,value=Efreq  if Efreq  is not None else 0.0 )
        dself.Ephase = depend_value(name="Ephase",value=Ephase if Ephase is not None else 0.0 )
        dself.Epeak  = depend_value(name="Epeak" ,value=Epeak  if Epeak  is not None else 0.0)
        dself.Esigma = depend_value(name="Esigma",value=Esigma if Esigma is not None else np.inf)
        dself.bec    = depend_array(name="bec"   ,value=bec    if bec    is not None else np.zeros(0))
        dself.cbec   = depend_array(name="cbec"   ,value=cbec)
        dself.cdip   = depend_array(name="cdip"   ,value=cdip)
        # dself.tacc   = depend_array(name="tacc"   ,value=tacc)

        # try : 
        #     self.eda = EDA(Eamp,Efreq,Ephase,Epeak,Esigma,cdip,cbec,bec)
        # except :
        #     print("WARNING: 'EDA' object not instantiated")

    def copy(self):
        return Ensemble(
            eens=self.eens,
            econs=0.0,
            temp=self.temp,
            pext=self.pext,
            stressext=dstrip(self.stressext).copy(),
            bcomponents=self.bcomp,
            bweights=dstrip(self.bweights).copy(),
            hweights=dstrip(self.hweights).copy(),
            time=self.time,
            Eamp=self.Eamp,
            Efreq=self.Efreq,
            Ephase=self.Ephase,
            Epeak=self.Epeak,
            Esigma=self.Esigma,
            bec=self.bec,
            cdip=self.cdip,
            # tacc=self.tacc,
            cBEC=self.cBEC,
            )

    def bind(
        self,
        beads,
        nm,
        cell,
        bforce,
        fflist,
        output_maker,
        enstype,
        elist=[],
        xlpot=[],
        xlkin=[],
    ):
        self.beads = beads
        self.cell = cell
        self.forces = bforce
        self.nm = nm
        dself = dd(self)
        self.output_maker = output_maker

        # this binds just the explicit bias forces
        self.bias.bind(
            self.beads,
            self.cell,
            self.bcomp,
            fflist,
            open_paths=nm.open_paths,
            output_maker=self.output_maker,
        )

        dself.econs = depend_value(name="econs", func=self.get_econs)
        # dependencies of the conserved quantity
        dself.econs.add_dependency(dd(self.nm).kin)
        dself.econs.add_dependency(dd(self.forces).pot)
        dself.econs.add_dependency(dd(self.bias).pot)
        dself.econs.add_dependency(dd(self.nm).vspring)
        dself.econs.add_dependency(dself.eens)

        # pipes the weights to the list of weight vectors
        i = 0
        for fc in self.bias.mforces:
            if fc.weight != 1:
                warning(
                    "The weight given to forces used in an ensemble bias are given a weight determined by bias_weight"
                )
            dpipe(dself.bweights, dd(fc).weight, i)
            i += 1

        # add Hamiltonian REM bias components
        if len(self.hweights) == 0:
            self.hweights = np.ones(len(self.forces.mforces))

        dself.hweights = depend_array(name="hweights", value=np.asarray(self.hweights))

        # we use ScaledForceComponents to replicate the physical forces without (hopefully) them being actually recomputed
        for ic in range(len(self.forces.mforces)):
            sfc = ScaledForceComponent(self.forces.mforces[ic], 1.0)
            self.bias.add_component(self.forces.mbeads[ic], self.forces.mrpc[ic], sfc)
            dd(sfc).scaling._func = lambda i=ic: self.hweights[i] - 1
            dd(sfc).scaling.add_dependency(dself.hweights)

        self._elist = []

        for e in elist:
            self.add_econs(e)

        dself.lpens = depend_value(
            name="lpens", func=self.get_lpens, dependencies=[dself.temp]
        )
        dself.lpens.add_dependency(dd(self.nm).kin)
        dself.lpens.add_dependency(dd(self.forces).pot)
        dself.lpens.add_dependency(dd(self.bias).pot)
        dself.lpens.add_dependency(dd(self.nm).vspring)

        # extended Lagrangian terms for the ensemble
        self._xlpot = []
        for p in xlpot:
            self.add_xlpot(p)

        self._xlkin = []
        for k in xlkin:
            self.add_xlkin(k)

        dself.eda.bind(self,enstype)

        
    def add_econs(self, e):
        self._elist.append(e)
        dd(self).econs.add_dependency(e)

    def add_xlpot(self, p):
        self._xlpot.append(p)
        dd(self).lpens.add_dependency(p)

    def add_xlkin(self, k):
        self._xlkin.append(k)
        dd(self).lpens.add_dependency(k)

    def get_econs(self):
        """Calculates the conserved energy quantity for constant energy
        ensembles.
        """

        eham = self.nm.vspring + self.nm.kin + self.forces.pot

        eham += self.bias.pot  # bias

        for e in self._elist:
            eham += e.get()

        return eham + self.eens

    def get_lpens(self):
        """Returns the ensemble probability (modulo the partition function)
        for the ensemble.
        """

        lpens = self.forces.pot + self.bias.pot + self.nm.kin + self.nm.vspring

        # inlcude terms associated with an extended Lagrangian integrator of some sort
        for p in self._xlpot:
            lpens += p.get()
        for k in self._xlkin:
            lpens += k.get()

        lpens *= -1.0 / (Constants.kb * self.temp * self.beads.nbeads)
        return lpens
        
class BEC(dobject):

    def __init__(self,cbec,bec):
        self.cbec = cbec
        dd(self).bec  = depend_array(name="bec",value=bec)
        pass

    def bind(self,eda,ensemble,enstype):

        self.enstype = enstype
        self.nbeads  = ensemble.beads.nbeads
        self.natoms  = ensemble.beads.natoms
        self.forces  = ensemble.forces # is this a weakref??
        # self.first   = True

        dself = dd(self)        
        if self.cbec: 
            dself.bec = depend_array(name="bec",\
                                     value=np.full((self.nbeads,3*self.natoms,3),np.nan),\
                                     # value=np.full((self.natoms,3,3),np.nan),\
                                     func=self._get_otf_BEC,\
                                     dependencies=[dd(eda).time,dd(ensemble.beads).q]) 
        elif self.enstype in EDA.integrators :
            temp = self._get_static_BEC() # reshape the BEC once and for all
            dself.bec = depend_array(name="bec",value=temp)
        else :
            # dself.bec = depend_array(name="bec",value=np.full((self.natoms,3,3),np.nan)) 
            dself.bec = depend_array(name="bec",value=np.full((self.nbeads,3*self.natoms,3),np.nan)) 

        # self.first = False
        pass

    def store(self,bec):
        super(BEC,self).store(bec)
        self.cbec.store(bec.cbec)

    def _get_otf_BEC(self,bead=None):
        """Return the BEC tensors (in cartesian coordinates), when computed on the fly (otf) by the driver"""
        # self._check_BEC()#skip=self.first)

        # print("using '_get_otf_BEC'")

        msg = "Error in '_get_otf_BEC'"

        # check that bead is a correct value
        if bead is not None:
            if bead < 0:
                raise ValueError("Error in '_get_otf_BEC': 'bead' is negative") 
            if bead >= self.nbeads :
                raise ValueError("Error in '_get_otf_BEC': 'bead' is greater than the number of beads") 
        else :
            if self.nbeads != 1 :
                raise ValueError("Error in '_get_otf_BEC': EDA integration has not implemented yet for 'nbeads' > 1")

        if self.cbec :
            if "BEC" not in self.forces.extras :
                raise ValueError(msg+": BEC tensors are not returned to i-PI (or at least not accessible in '_get_otf_BEC').") 
        else :
            raise ValueError(msg+": you should not get into this functon if 'cbec' is False.") 
        
        BEC = np.full((self.nbeads,3*self.natoms,3),np.nan)
        for n in range(self.nbeads):
            bec = np.asarray(self.forces.extras["BEC"][n])

            if bec.shape[0] != 3 * self.natoms :
                raise ValueError(msg+": number of BEC tensors is not equal to the number fo atoms x 3.")
            if bec.shape[1] != 3 :
                raise ValueError(msg+": BEC tensors with wrong shape. They should have 3 components.")
            
            BEC[n,:,:] = np.copy(bec)
            
            # Na = self.natoms
            # bec = bec.reshape((Na,3,3))
            # Axis of bec :
            #   1st: atoms index (0,1,2...)
            #   2nd: atom coordinate (x,y,z)
            #   3rd: dipole direction (x,y,z)
        return BEC

    def _get_static_BEC(self):
        """Return the BEC tensors (in cartesian coordinates).
        The BEC tensor are stored in a compact form.
        This method trasform the BEC tensors into another data structure, suitable for computation.
        A lambda function is also returned to perform fast matrix multiplication.
        """

        raise ValueError("This function has to be re-written")
        # self.first = True 

        N = len(self.bec)      # lenght of the BEC array
        Na = self.natoms # number of atoms

        if N == Na:     # scalar BEC
            Z = np.zeros((Na,3,3))
            for i in range(Na):
                for j in range(3):
                    Z[i,j,j] = self.bec[i]
                Z[i,:,:] = Z[i,:,:]
            return Z #self._lv2cart(Z)
            #lambda a,b : a*b # element-wise (matrix) multplication (only the diagonal elements have been allocated)

        elif N == 3*Na: # diagonal BEC
            Z = np.zeros((Na,3,3))
            temp = self.bec.reshape((Na,3))
            for i in range(Na):
                for j in range(3):
                    Z[i,j,j] = temp[i,j]
                Z[i,:,:] = Z[i,:,:]
            return # self._lv2cart(Z)
            #lambda a,b : a*b # element-wise (matrix) multplication (only the diagonal elements have been allocated)
        
        elif N == 9*Na: # all-components BEC
            Z = np.zeros((Na,3,3))
            temp = self.bec.reshape((Na,3,3))
            for i in range(Na):
                Z[i,:,:] = temp[i,:,:]
            return Z # self._lv2cart(Z) # rows-by-columns (matrix) multplication (all the elements have been allocated)
        # elif N == 0 and empty_is_okay:
        #     return np.zeros((Na,3,3))
        else :
            raise ValueError("BEC tensor with wrong size!")

    # def _check_BEC(self):#,skip=False):
    #     """Check that the BEC tensors are correctly formatted."""

    #     if self.nbeads != 1 :
    #         raise ValueError("Error in '_check_BEC': EDA integration has not implemented yet for 'nbeads' > 1")
        
    #     msg = "Error in '_check_BEC'"

    #     if self.cbec :
    #         if "BEC" not in self.forces.extras :
    #             raise ValueError(msg+": BEC tensors are not returned to i-PI (or at least not accessible in '_check_BEC').") 
    #     else :
    #         return True

    #     bec = np.asarray(self.forces.extras["BEC"][0])

    #     if bec.shape[0] != self.natoms :
    #         raise ValueError(msg+": number of BEC tensors is not equal to the number fo atoms.")
    #     if bec.shape[1] != 9 :
    #         raise ValueError(msg+": BEC tensors with wrong shape. They should have 9 components.")

    #     return True
    
    #     # if len(self.forces.extras["BEC"]) != Nb:
    #     #     raise ValueError(msg+": wrong number of bead for the BEC tensors.")
        
    #     # check whether the BEC tensors have the correct shape 
    #     becs = self.forces.extras["BEC"]
    #     for i in range(Nb):
    #         Na = len(becs[i])
    #         if Na != self.natoms:
    #             raise ValueError(msg+": number of BEC tensors is not equal to the number fo atoms.")
    #         bec = becs[i]
    #         for j in range(Na):
    #             if len(bec[j]) != 9 :
    #                 raise ValueError(msg+": BEC tensors with wrong shape. They should have 9 components.")
    #     return True

class Dipole(dobject):

    def __init__(self,cdip):
        self.cdip = cdip
        pass

    def bind(self,eda,ensemble):        
        dself = dd(self)
        dself.nbeads = depend_value(name="nbeads",value=ensemble.beads.nbeads ) 

        self.forces = ensemble.forces # is this a weakref??

        val = np.full(self.nbeads,np.zeros(3,dtype=float)) if self.nbeads > 1 else np.zeros(3,dtype=float)
        dself._dipole_ = depend_array(name="_dipole_", func=lambda:self._get_dipole(bead=0),value=val,dependencies=[dd(eda).time,dd(ensemble.beads).q])

        pass

    def store(self,dipole):
        super(Dipole, self).store(dipole)
        self.cdip.store(dipole.cdip)
        pass

    def _get_dipole(self,bead=None):
        """Return the electric dipole of all the beads as a list of np.array"""
        # self._check_dipole()

        # check that bead is a correct value
        # N = self.beads.nbeads
        if bead is not None:
            if bead < 0:
                raise ValueError("Error in '_get_dipole': 'beads' is negative.") 
            if bead >= self.nbeads :
                raise ValueError("Error in '_get_dipole': 'beads' is greater than the number of beads.") 
            if bead > 1 :
                raise ValueError("The case with 'beads' != 0 has not been implemeted yet") 
            
        if not self.cdip:
            return np.asarray([0,0,0])
        else :
            if "dipole" in self.forces.extras :
                dipole = np.asarray(self.forces.extras["dipole"]).flatten()
                if len(dipole) != 3 :
                    print("dipole:",dipole)
                    raise ValueError("'dipole' has not length 3")
                return dipole
                # dipole = [ self.forces.extras["dipole"][i] for i in range(self.nbeads)]
                # return dipole[0] if bead is None else dipole[bead] 

            elif "raw" not in self.forces.extras :
                raise ValueError("'raw' has to be in 'forces.extras'")

            elif np.all( [ "Total dipole moment" in s for s in self.forces.extras["raw"] ] )  :
                raw = [ self.forces.extras["raw"][i] for i in range(self.nbeads)]
                raw = raw[0] if bead is None else raw[bead] 
                factor = 1.
                if "[eAng]" in raw : 
                    factor = UnitMap["length"]["angstrom"]

                dipole = np.full(3,np.nan)
                pattern = "[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|\b[-+]?\d+\b"
                matches = re.findall(pattern, raw)
                if len(matches) != 3 :
                    raise ValueError("wrong number of extracted values from the extra string: they should be 3.")
                else :
                    for n in range(3):
                        dipole[n] = float(matches[n])
                return float(factor) * dipole
            else :
                raise ValueError("Error in '_get_dipole': can not extract dipole from the extra string.") 
               
        
    # def _check_dipole(self):
    #     """Check that the electric dipole is correctly formatted."""

    #     print("ELIA ->           type(self.forces.extras): ",type(self.forces.extras))
    #     print("ELIA ->                 self.forces.extras: ",self.forces.extras)
    #     print("ELIA ->          self.forces.extras.keys(): ",self.forces.extras.keys())
    #     print("ELIA -> self.forces.extras[\'raw\']: ",self.forces.extras["raw"])
        
    #     msg = "Error in '_check_dipole'"
    #     error = ValueError(msg+": the dipole is not returned to i-PI (or at least not accessible in '_check_dipole').")

    #     if self.cdip :

    #         if "dipole" in self.forces.extras :
    #             arr = np.asrray(self.forces.extras["dipole"])
    #             print("ELIA -> dipole.shape :",arr.shape)
    #             if len(self.forces.extras["dipole"]) != self.nbeads:
    #                 raise ValueError(msg+": wrong number of bead for the dipole.")

    #         else :
    #             if "raw" not in self.forces.extras :
    #                 raise error
    #             else :
    #                 if len(self.forces.extras["raw"]) != self.nbeads:
    #                     raise ValueError(msg+": wrong number of bead for extra strings.")
                    
    #                 if not np.all( [ "Total dipole moment" in s for s in self.forces.extras["raw"] ] ) :
    #                     raise error
    #             return True


    #     else :
    #         return True       


class ElectricField(dobject):

    def __init__(self,Eamp,Efreq,Ephase,Epeak,Esigma):

        dself = dd(self)
        dself.Eamp   = depend_array(name="Eamp"  ,value=Eamp   if Eamp   is not None else np.zeros(3))
        dself.Efreq  = depend_value(name="Efreq" ,value=Efreq  if Efreq  is not None else 0.0 )
        dself.Ephase = depend_value(name="Ephase",value=Ephase if Ephase is not None else 0.0 )
        dself.Epeak  = depend_value(name="Epeak" ,value=Epeak  if Epeak  is not None else 0.0)
        dself.Esigma = depend_value(name="Esigma",value=Esigma if Esigma is not None else np.inf)     

        # these will be overwritten
        dself.Eenvelope      = depend_value(name="Eenvelope" ,value=1.0,func=self._get_Eenvelope)
        dself.TderEenvelope  = depend_value(name="TderEenvelope" ,value=0.0,func=self._get_TderEenvelope)  
        dself.Efield         = depend_array(name="Efield",    value=np.zeros(3, float),func=self._get_Efield)
        dself.TderEfield     = depend_array(name="TderEfield",value=np.zeros(3, float),func=self._get_TderEfield)

    def bind(self,eda,enstype):

        self.enstype = enstype
        dself = dd(self)
        dself.cptime = depend_value(name="cptime",value=0.0)
        dpipe(dfrom=dd(eda).cptime,dto=dd(self).cptime)

        # same dependencies for Eenvelope and its time derivative
        dep = [dself.cptime,dself.Epeak,dself.Esigma] 
        dself.Eenvelope      = depend_value(name="Eenvelope" ,value=1.0,func=self._get_Eenvelope,dependencies=dep)
        dself.TderEenvelope  = depend_value(name="TderEenvelope" ,value=0.0,func=self._get_TderEenvelope,dependencies=dep)

         
        if enstype in EDA.integrators:
            # with dependencies
            dep = [dself.cptime,dself.Eamp,dself.Efreq,dself.Ephase,dself.Eenvelope]
            dself.Efield     = depend_array(name="Efield",    value=np.zeros(3, float),func=self._get_Efield,    dependencies=dep)
            dself.TderEfield = depend_array(name="TderEfield",value=np.zeros(3, float),func=self._get_TderEfield,dependencies=dep)
        else :
            # no dependencies
            dself.Efield     = depend_array(name="Efield",    value=np.zeros(3, float),func=lambda time=None : np.zeros(3, float))
            dself.TderEfield = depend_array(name="TderEfield",value=np.zeros(3, float),func=lambda time=None : np.zeros(3, float))
        
        pass

    def store(self,ef):
        super(ElectricField, self).store(ef)
        self.Eamp.store(ef.Eamp)
        self.Efreq.store(ef.Efreq)
        self.Ephase.store(ef.Ephase)
        self.Epeak.store(ef.Epeak)
        self.Esigma.store(ef.Esigma)
        pass

    # def __get__(self):
    #     return self.Efield
    
    # def __call__(self,time=None):
    #     # bypass self.Efield
    #     return self._get_Efield(time)
    
    def _get_Efield(self,time=None):
        """Get the value of the external electric field (cartesian axes)"""
        if time is None :
            raise ValueError("Hey man! Don't you think it's better to specify the time you want to evaluate the electric field?")
        if hasattr(time, "__len__"):
            return np.outer ( self._get_Ecos(time) * dd(self).Eenvelope(time) , self.Eamp )
        else :
            return self._get_Ecos(time) * dd(self).Eenvelope(time) * self.Eamp
         
    
    def _Eenvelope_is_on(self):
        return self.Epeak > 0.0 and self.Esigma != np.inf
    
    def _get_Eenvelope(self,time=None):
        """Get the gaussian envelope function of the external electric field"""
        # https://en.wikipedia.org/wiki/Normal_distribution
        if self._Eenvelope_is_on() :
            x = self.cptime if time is None else time # indipendent variable
            u = self.Epeak  # mean value
            s = self.Esigma # standard deviation
            return np.exp( - 0.5 * ((x-u)/s)**2 ) # the returned maximum value is 1, when x = u
        else :
            return 1.0
        
    def _get_TderEenvelope(self,time=None):
        """Get the time derivative of Eenvelope"""
        # https://en.wikipedia.org/wiki/Normal_distribution
        if self._Eenvelope_is_on() :
            x = self.cptime if time is None else time # indipendent variable
            u = self.Epeak  # mean value
            s = self.Esigma # standard deviation
            return - dd(self).Eenvelope(time) * ( x - u ) / s**2
        else :
            return 0.0
            #raise warning("Time derivative of Eenvelope should be evaluate only when Eenvelope is non vanishing")

    def _get_TderEfield(self,time=None):
        """Get the time derivative of the external electric field (cartesian axes)"""
        #if self.Eenvelope is not None : 
        #  chain rule
        if time is None :
            raise ValueError("Hey man! Don't you think it's better to specify the time you want to evaluate the (time derivative of the) electric field?")
        return self.Eamp * ( self._get_Ecos(time) * self._get_TderEenvelope(time) + self._get_TderEcos(time) * dd(self).Eenvelope(time) ) 
        #else :
        #    return self.Eamp * self._get_TderEcos()

    def _get_Ecos(self,time=None):
        """Get the sinusoidal part of the external electric field"""
        t = self.cptime if time is None else time
        return np.cos( self.Efreq * t + self.Ephase)
    
    def _get_TderEcos(self,time=None):
        """Get the time derivative sinusoidal part of the external electric field"""
        t = self.cptime if time is None else time
        return - self.Efreq * np.sin( self.Efreq * t + self.Ephase)
    
class EDA(dobject):

    integrators = ["eda-nve","eda-nvt"]

    def __init__(self,Eamp,Efreq,Ephase,Epeak,Esigma,cdip,cbec,bec,**kwargv):
        super(EDA,self).__init__(**kwargv)
        self.Electric_Field = ElectricField(Eamp,Efreq,Ephase,Epeak,Esigma)
        self.Dipole         = Dipole(cdip)
        self.Born_Charges   = BEC(cbec,bec)
        # self.tacc           = depend_value(name="tacc" ,value=tacc)
        pass

    def bind(self,ensemble,enstype):

        self.enstype = enstype
        dself = dd(self)

        dself.econs  = depend_value(name="econs",value=0.0)
        dself.time   = depend_value(name="time",value=0.0)
        dself.cptime = depend_value(name="cptime",value=0.0)

        dpipe(dfrom=dd(ensemble).econs, dto=dself.econs)
        dpipe(dfrom=dd(ensemble).time,  dto=dself.time)

        self.Electric_Field.bind(self,enstype)
        self.Dipole.bind(self,ensemble)
        self.Born_Charges.bind(self,ensemble,enstype)  
        
        # for easier access
        dself.Efield   = depend_array(name="Efield",value=np.full(dd(self.Electric_Field).Efield.shape,np.nan),func=lambda time=None: dd(self.Electric_Field).Efield(time))
        dself.bec      = depend_array(name="bec",   value=np.full(dd(self.Born_Charges).bec.shape,     np.nan))
        dself.dipole   = depend_array(name="dipole",value=np.full(dd(self.Dipole)._dipole_.shape,np.nan))

        dpipe(dfrom=dd(self.Born_Charges).bec,     dto=dself.bec)
        dpipe(dfrom=dd(self.Dipole)._dipole_,dto=dself.dipole)  

        # experimental
        dpipe(dfrom=dself.time,dto=dself.cptime)

        dep = [dself.dipole,dself.Efield,dself.time]
        # dself.EDAenergy     = depend_value(name="EDAenergy",     func=self._get_EDAenergy,    value=0.0,dependencies=dep)
        # dself.TderEDAenergy = depend_value(name="TderEDAenergy", func=self._get_TderEDAenergy,value=0.0,dependencies=dep)
        # dself.Eenthalpy     = depend_value(name="Eenthalpy",     func=self._get_Eenthalpy,    value=0.0,dependencies=dep)

        # dself.tacc
        # dself.Tconserved = depend_value(name="Tconserved", func=self._get_Tconserved,value=0.0,dependencies=[dself.Eenthalpy,dself.time,dself.cptime])

        pass

    def store(self,eda):
        super(EDA, self).store(eda)
        self.Electric_Field.store(eda.Electric_Field)
        self.Dipole.store(eda.Dipole)
        self.Born_Charges.store(eda.bec)
        # self.tacc.store(eda.tacc)
        pass        

    # def _get_EDAenergy(self,time=None):
    #     """EDA contribution to the enthalpy"""
    #     return float( np.dot( self.dipole , dd(self.Electric_Field).Efield(time) ))
    
    # def _get_TderEDAenergy(self,time=None): # ES: to be modified
    #     """Time derivative of EDAenergy"""
    #     #self._check_time(msg="calling")
    #     return float( np.dot( self.dipole , dd(self.Electric_Field).TderEfield(time) ))

    # def _get_Eenthalpy(self,time=None):
    #     """Electric enthalpy"""
    #     return self.econs - dd(self).EDAenergy(time)
              
    # def _get_Tconserved(self,time=None):
    #     """Conserved quantity for time-dependent systems"""
    #     # tacc is added and not subtracted because it is defined using EDAenergy
    #     # but EDAenergy is subtracted the energy
    #     # # so we have two minus signs 
    #     self._check_time(msg="calling",time=time)
    #     return dd(self).Eenthalpy(time) + self.tacc

    def _check_time(self,msg="coding",time=None):
        """Check that self.cptime is equal to self.ensemble.time.
        Pay attention that this is not always true all over the simulation!
        These variable have to be equal only before and after the Integration procedure.
        In fact, this method is called only in Dynamics.step, after self.integrator.step(step).
        The two variable are also forces to be equal before the INtegration procedure at each step.

        This method should always return True, but perhaps future code changes could "break" this.
        Better to be sure that everythin is fine :) """

        coding = "Error in Ensemble._check_time: the 'continous' time 'Ensemble.cptime' does not match"+\
                "Ensemble.time (up to a threshold).\nThis seems to be a coding error, not due to wrong input parameters."+\
                "\nRead the description of the function in file ipi/engine/motion/dynamics.py."+\
                "\nAnd then, if you still have problem, you can write me an email to stocco@fhi-berlin.mpg.de.\nBye :)"
        
        calling = "Error in Ensemble._check_time: the 'continous' time 'Ensemble.cptime' does not match"+\
                "Ensemble.time (up to a threshold).\nThis seems to be a coding error, not due to wrong input parameters."+\
                "\nIt seems that you are trying to evaluate the time derivative of some variables in the wrong moment."+\
                "\nRead the description of the function in file ipi/engine/motion/dynamics.py."+\
                "\nAnd then, if you still have problem, you can write me an email to stocco@fhi-berlin.mpg.de.\nBye :)"
        
        messages = {"coding":coding,"calling":calling}
        
        # if the specified msg does not exists, the code will tell you and it will print the 'coding' error message
        if msg not in messages.keys():
            messages[msg] = "'{:s}' is not a valid message key for _check_time.\n".format(msg) + messages["coding"]

        thr_time_comparison = 0.1
        if abs(self.cptime - self.time) > thr_time_comparison:
            raise ValueError(messages[msg])
        if time is not None and abs(time - self.time) > thr_time_comparison:
            raise ValueError(messages[msg])
        return True
