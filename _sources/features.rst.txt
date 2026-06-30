Core features
=============

i-PI includes a large number of advanced molecular dynamics features,
with an obvious focus on path integral molecular dynamics, but also
several methods for sampling classical trajectories. 
This is an (incomplete) list of some of the main features in alphabetical order.

If you implement a major new feature, this is the place to briefly outline 
what it does and what are the papers to cite. Please follow as closely as 
possible the template.

Bosonic and Fermionic Path Integral Molecular Dynamics
------------------------------------------------------

PIMD simulations of bosonic particles with a polynomial scaling
algorithm. Supports mixtures of distinguishable and bosonic particles.
Fermionic statistics can be obtained by a reweigthing procedure by post
processing the simulation (see ref. 2 below).

| **Main contributors:** Yotam Feldman, Barak Hirshberg
| **Theory and Implementation:**
| Y.M.Y. Feldman and B. Hirshberg *“Quadratic scaling bosonic path
  integral molecular dynamics simulations”*, J. Chem. Phys. 159, 154107 (2023) DOI: 
  `10.1063/5.0173749 <https://doi.org/10.1063/5.0173749>`__
  — BibTeX:
  `fetch <https://pubs.aip.org/Citation/Download?resourceId=2917443&resourceType=3&citationFormat=2>`__
| **Theory:**
| B. Hirshberg, V. Rizzi and M. Parrinello, *“Path integral molecular
  dynamics for bosons”*, Proc. Natl. Acad. Sci. U.S.A. 116 (43)
  21445-21449 (2019) DOI:
  `10.1073/pnas.1913365116 <https://doi.org/10.1073/pnas.1913365116>`__
  — BibTeX:
  `fetch <https://www.doi2bib.org/bib/10.1073/pnas.1913365116>`__
| B. Hirshberg, M. Invernizzi and M. Parrinello, *“Path integral
  molecular dynamics for fermions: Alleviating the sign problem with the
  Bogoliubov inequality”*, J. Chem. Phys. 152, 171102 (2020) DOI:
  `10.1063/5.0008720 <https://doi.org/10.1063/5.0008720>`__ — BibTeX:
  `fetch <https://www.doi2bib.org/bib/10.1063/5.0008720>`__
| C. W. Myung, B. Hirshberg and M. Parrinello, *“Prediction of a
  Supersolid Phase in High-Pressure Deuterium”*, Phys. Rev. Lett., 128
  045301 (2022) DOI:
  `10.1103/PhysRevLett.128.045301 <https://doi.org/10.1103/PhysRevLett.128.045301>`__
  — BibTeX:
  `fetch <https://www.doi2bib.org/bib/10.1103/PhysRevLett.128.045301>`__


Cavity Molecular Dynamics for Polaritonics
------------------------------------------

This initial implementation provides an efficient cavity molecular
dynamics (CavMD) scheme for simulating strong light-matter interactions
between molecules and an optical cavity mode, particularly in the
vibrational strong coupling regime. At present, the nuclear partial
charges are assumed to be fixed during the simulation. CavMD is
implemented with a new force evaluator: *ffcavphsocket*, which is
operated similarly to the original *ffsocket* evaluator, but with
additional parameters for controlling cavity photons. Hence, with this
implementation, users can study different aspects of vibrational strong
coupling with many sophisticated methods supported in i-pi.

| **Main contributors:** Tao E. Li
| **Implementation and theory:**
| T. E. Li, J. E. Subotnik, and A. Nitzan, *“Cavity molecular dynamics
  simulations of liquid water under vibrational ultrastrong coupling”*,
  Proc. Natl. Acad. Sci. 117(31), 18324–18331. (2020)
| DOI:
  `10.1073/pnas.2009272117 <http://dx.doi.org/10.1073/pnas.2009272117>`__
  — BibTeX:
  `fetch <https://www.doi2bib.org/bib/10.1073/pnas.2009272117>`__
| T. E. Li, A. Nitzan, S. Hammes-Schiffer, and J. E. Subotnik, *“Quantum
  simulations of vibrational strong coupling via path integrals”*, J.
  Phys. Chem. Lett. 13(17), 3890–3895. (2022)
| DOI:
  `10.1021/acs.jpclett.2c00613 <http://dx.doi.org/10.1021/acs.jpclett.2c00613>`__
  — BibTeX:
  `fetch <https://www.doi2bib.org/bib/10.1021/acs.jpclett.2c00613>`__

Committee models
----------------

Uncertainty estimation for machine-learning potentials based on
committee models. Multiple potentials are used simultaneously (they
should have been fitted with randomized training sets). The mean is used
to drive the system, the spread is used as an estimate of the
uncertainty. This uncertainty can be used to select high-error
structures for active learning, to estimate the propagation of the error
to thermodynamic averages, or to build a weighted baseline model that
falls back to a safe baseline potential when the ML models fail.

| **Main contributors:** Giulio Imbalzano, Venkat Kapil, Yongbin Zhuang,
  Federico Grasselli, Michele Ceriotti
| **Implementation and theory:**
| G. Imbalzano, Y. Zhuang, V. Kapil, K. Rossi, E. A. Engel, F.
  Grasselli, and M. Ceriotti, *“Uncertainty estimation for molecular
  dynamics and sampling”*, J. Chem. Phys. 154(7), 074102 (2021)
| DOI: `10.1063/5.0036522 <http://dx.doi.org/10.1063/5.0036522>`__ — BibTeX:
  `fetch <http://www.doi2bib.org/#/doi/10.1063/5.0036522>`__
| F. Musil, M. J. Willatt, M. A. Langovoy, and M. Ceriotti, *“Fast and
  Accurate Uncertainty Estimation in Chemical Machine Learning”*,
  Journal of Chemical Theory and Computation 15(2), 906–915 (2019)
| DOI:
  `10.1021/acs.jctc.8b00959 <http://dx.doi.org/10.1021/acs.jctc.8b00959>`__ —
  BibTeX:
  `fetch <http://www.doi2bib.org/#/doi/10.1021/acs.jctc.8b00959>`__

Direct Estimators for Isotope Fractionation
-------------------------------------------

A direct estimator to evaluate the isotope fractionation ratios using a
single operation (and a single keyword in the input file), without the
need for a thermodynamic integration with respect to the mass of the
isotope.

| **Main contributors:** Bingqing Cheng, Michele Ceriotti
| **Implementation and Theory:**
| B.Cheng, M.Ceriotti, “Direct path integral estimators for isotope
  fractionation ratios.” The Journal of chemical physics 141, 244112
  (2015)
| DOI: `10.1063/1.4904293 <http://dx.doi.org/10.1063/1.4904293>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.4904293>`__

Fast-Forward Langevin Thermostat
--------------------------------

This is a modified form of Langevin dynamics in which sluggish
high-friction behaviour is corrected for by flipping a particle’s
momentum when the action of the thermostat causes it to change
direction.

| **Main contributors:** Mahdi Hijazi, David Wilkins, Michele Ceriotti
| **Implementation and Theory:**
| M. Hijazi, D. M. Wilkins, *“Fast-forward Langevin dynamics with
  momentum flips”*, J. Chem. Phys. 148, 184109 (2018)
| DOI: `10.1063/1.5029833 <http://dx.doi.org/10.1063/1.5029833>`__ — BibTeX:
  `fetch <https://www.doi2bib.org/bib/10.1063%2F1.5029833>`__


Finite-differences Suzuki-Chin PIMD
-----------------------------------

Suzuki-Chin PIMD gives better convergence w.r.t. the number of imaginary
time slices as compared to the standard Trotter scheme. The
implementation uses a symplectic and time-reversible finite-difference
algorithm to compute high order corrections to traditional PIMD for any
empirical or ab initio forcefield.

| **Main contributors:** Venkat Kapil, Michele Ceriotti
| **Implementation and Theory:**
| V.Kapil, J.Behler, M.Ceriotti *“High order path interals made easy”*,
  J. Chem. Phys. 145, 234103 (2016)
| DOI: `10.1063/1.4971438 <http://dx.doi.org/10.1063/1.4971438>`__ —
  BIBTEX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.4971438>`__
| **Theory:**
| S.Jang, S.Jang, G.A.Voth *“Applications of higher order composite
  factorization schemes in imaginary time path integral simulations”*,
  J. Chem. Phys. 115, 7832 (2001)
| DOI: `10.1063/1.1410117 <http://dx.doi.org/10.1063/1.1410117>`__ —
  BIBTEX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.1410117>`__
| S.A.Chin *“Symplectic integrators from composite operator
  factorizations”*, Phys. Lett. A 226, 344 (1997)
| DOI:
  `10.1016/s0375-9601(97)00003-0 <http://dx.doi.org/10.1016/s0375-9601(97)00003-0>`__
  — BIBTEX:
  `fetch <http://www.doi2bib.org/bib/10.1016%2Fs0375-9601(97)00003-0>`__
| M.Suzuki *“Hybrid exponential product formulas for unbounded operators
  with possible applications to Monte Carlo simulations”*, Phys. Lett. A
  201, 425 (1995)
| DOI:
  `10.1016/0375-9601(95)00266-6 <http://dx.doi.org/10.1016/0375-9601(95)00266-6>`__
  — BIBTEX:
  `fetch <http://www.doi2bib.org/bib/10.1016%2F0375-9601(95)00266-6>`__

Finite-differences Vibrational Analysis
---------------------------------------

Harmonic vibrations through finite differences for simple evaluation of
the harmonic Hessian.

| **Main contributors:** Kapil, Bienvenue
| **Implementation:**
| M. Rossi, P. Gasparotto, M. Ceriotti, *“Anharmonic and Quantum
  Fluctuations in Molecular Crystals: A First-Principles Study of the
  Stability of Paracetamol”*, Phs. Rev. Lett. 117, 115702 (2016)
| DOI:
  `10.1103/PhysRevLett.117.115702 <http://dx.doi.org/10.1103/PhysRevLett.117.115702>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/#/doi/10.1103/PhysRevLett.117.115702>`__

Free-energy Perturbation Estimators for Isotope Fractionation
-------------------------------------------------------------

Computing isotope fractionation using the thermodynamic integration
method requires evaluating the quantum kinetic energy of several systems
containing atoms that have different fictitious masses between the
physical masses of two isotopes, meaning that a number of PIMD
simulations have to be performed. With the help of re-weighting, one has
the option of running just one set of simulation with a certain
fictitious mass, and obtain the quantum kinetic energy for systems with
other masses.

| **Main contributors:** Michele Ceriotti, Thomas Markland
| **Theory and implementation:**
| Michele Ceriotti, Thomas E. Markland, “Efficient methods and practical
  guidelines for simulating isotope effects.” The Journal of chemical
  physics 138(1), 014112 (2013).
| DOI: `10.1063/1.4772676 <http://dx.doi.org/10.1063/1.4772676>`__ — BibTeX:
  `fetch <http://www.doi2bib.org/#/doi/10.1063/1.4772676>`__

Generalized Langevin Equation Thermostats
-----------------------------------------

The Generalized Langevin Equation provides a very flexible framework to
manipulate the dynamics of a classical system, improving sampling
efficiency and obtaining quasi-equilibrium ensembles that mimic quantum
fluctuations. Parameters for the different modes of operation can be
obtained from the `GLE4MD
website <http://gle4md.org/index.html?page=matrix>`__.

| **Main contributors:** Michele Ceriotti
| **Implementation:**
| M. Ceriotti, G. Bussi, M. Parrinello, *“M. Colored-Noise Thermostats à
  la Carte”*, J. Chem. Theory Comput. 6, 1170–1180 (2010)
| DOI: `10.1021/ct900563s <http://dx.doi.org/10.1021/ct900563s>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1021%2Fct900563s>`__
| **Theory:**
| *Optimal Sampling Efficiency* — M. Ceriotti, G. Bussi, and M.
  Parrinello, *“Langevin Equation with Colored Noise for
  Constant-Temperature Molecular Dynamics Simulations”*, Phys.
  Rev. Lett. 102, 20601 (2009)
| DOI:
  `10.1103/PhysRevLett.102.020601 <http://dx.doi.org/10.1103/PhysRevLett.102.020601>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1103%2FPhysRevLett.102.020601>`__
| *Quantum Thermostat* — M. Ceriotti, G. Bussi, and M. Parrinello,
  *“Nuclear Quantum Effects in Solids Using a Colored-Noise
  Thermostat”*, Phys. Rev. Lett. 103, 30603 (2009)
| DOI:
  `10.1103/PhysRevLett.103.030603 <http://dx.doi.org/10.1103/PhysRevLett.103.030603>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1103%2FPhysRevLett.103.030603>`__
| *Delta Thermostat* — M. Ceriotti and M. Parrinello, *“The
  δ-Thermostat: Selective Normal-Modes Excitation by Colored-Noise
  Langevin Dynamics”*, Procedia Comput. Sci. 1, 1607 (2010)
| DOI:
  `10.1016/j.procs.2010.04.180 <http://dx.doi.org/10.1016/j.procs.2010.04.180>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1016%2Fj.procs.2010.04.180>`__
| *MTS Thermostat* — J. A. Morrone, T. E. Markland, M. Ceriotti, and B.
  J. Berne, *“Efficient Multiple Time Scale Molecular Dynamics: Using
  Colored Noise Thermostats to Stabilize Resonances”*, J. Chem. Phys.
  134, 14103 (2011)
| DOI: `10.1063/1.3518369 <http://dx.doi.org/10.1063/1.3518369>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.3518369>`__
| *“Hot-spot”* — R. Dettori, M. Ceriotti, J. Hunger, C. Melis, L.
  Colombo, and D. Donadio, *“Simulating Energy Relaxation in Pump-Probe
  Vibrational Spectroscopy of Hydrogen-Bonded Liquids”*, J. Chem. Theory
  Comput. (2017)
| DOI:
  `10.1021/acs.jctc.6b01108 <http://dx.doi.org/10.1063/10.1021/acs.jctc.6b01108>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1021%2Facs.jctc.6b01108>`__

Geometry Optimization
---------------------

Several standard algorithms for geometry optimization have been
implemented to give the convenience of static calculations that are
fully compatible with (PI)MD and other advanced sampling techniques.

| **Main contributors:** Benjamin Helfrecht, Sophie Mutzel, Riccardo
  Petraglia, Yair Litman, Mariana Rossi
| **Implementation:**
| M. Rossi, P. Gasparotto, and M. Ceriotti, *“Anharmonic and Quantum
  Fluctuations in Molecular Crystals: A First-Principles Study of the
  Stability of Paracetamol”*, Phys. Rev. Lett. 117, 115702 (2016)
| DOI:
  `10.1103/PhysRevLett.117.115702 <http://dx.doi.org/10.1103/PhysRevLett.117.115702>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1103%2FPhysRevLett.117.115702>`__
| **Theory:**
| W. H. Press, *“Numerical Recipes: The Art of Scientific Computing”*,
  (Cambridge University Press, 2007)

Langevin Sampling for Noisy or Dissipative Forces
-------------------------------------------------

A modified Langevin thermostat that allows for constant-temperature
dynamics with noisy or dissipative forces by applying additional damping
or noise for compensation. The implementation contains a method to
adjust the amount of compensation automatically.

| **Main contributors:** Jan Kessler, Thomas D. Kühne
| **Theory:**
| T. D. Kühne, M. Krack, F. R. Mohamed, M. Parrinello, *“Efficient and
  Accurate Car-Parrinello-like Approach to Born-Oppenheimer Molecular
  Dynamics”*, Phys. Rev. Lett. 98, 066401 (2007)
| DOI:
  `10.1103/PhysRevLett.98.066401 <http://dx.doi.org/10.1103/PhysRevLett.98.066401>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1103%2FPhysRevLett.98.066401>`__
| F. R. Krajewski, M. Parrinello, *“Linear scaling electronic structure
  calculations and accurate statistical mechanics sampling with noisy
  forces”*, Phys. Rev. B 73, 041105 (2006)
| DOI:
  `10.1103/PhysRevB.73.041105 <http://dx.doi.org/10.1103/PhysRevB.73.041105>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1103%2FPhysRevB.73.041105>`__
| Y. Luo, A. Zen, S. Sorella, *“Ab initio molecular dynamics with noisy
  forces: Validating the quantum Monte Carlo approach with benchmark
  calculations of molecular vibrational properties”*, J. Chem. Phys.
  141, 194112 (2014)
| DOI: `10.1063/1.4901430 <http://dx.doi.org/10.1063/1.4901430>`__ — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1063%2F1.4901430>`__

Multiple Time Step integrators
------------------------------

A multiple time step integration scheme allows for integration of
different components of forces with different time steps. It becomes
advantageous when the total force can be decomposed into a slowly
varying expensive part and a rapidly varying cheap part. A larger time
step can be used to integrate the former, there by reducing the number
of expensive computations.

| **Main contributors:** Venkat Kapil
| **Implementation:**
| V.Kapil, J.VandeVondele, M.Ceriotti *“Accurate molecular dynamics and
  nuclear quantum effects at low cost by multiple steps in real and
  imaginary time: using density functional theory to accelerate
  wavefunction methods”*, J. Chem. Phys. 144, 054111 (2016)
| DOI: `10.1063/1.4941091 <http://dx.doi.org/10.1063/1.4941091>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.4941091>`__
| **Theory:**
| M.Tuckerman, B.J.Berne *“Reversible multiple time scale molecular
  dynamics”*, J. Chem. Phys. 97, 1990 (1992)
| DOI: `10.1063/1.463137 <http://dx.doi.org/10.1063/1.463137>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.463137>`__

Open Path Integrals
-------------------

Open path integrals and momentum distribution estimators for the
computation of the particle momentum distribution including quantum
fluctuations of nuclei.

| **Main contributors:** Kapil, Cuzzocrea, Ceriotti
| **Implementation and Theory:**
| V. Kapil, A. Cuzzocrea, M. Ceriotti, *“Anisotropy of the Proton
  Momentum Distribution in Water”*, J. Phys. Chem. B 122, 6048-6054
  (2018)
| DOI:
  `10.1021/acs.jpcb.8b03896 <http://dx.doi.org/10.1021/acs.jpcb.8b03896>`__ —
  BibTeX:
  `fetch <http://www.doi2bib.org/#/doi/10.1021/acs.jpcb.8b03896>`__
| **Theory:**
| J. A. Morrone, R. Car, *“Nuclear Quantum Effects in Water”*, Phys.
  Rev. Lett. 101, 017801 (2008)
| DOI:
  `10.1103/PhysRevLett.101.017801 <http://dx.doi.org/10.1103/PhysRevLett.101.017801>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/#/doi/10.1103/PhysRevLett.101.017801>`__

Path Integral GLEs
------------------

Generalized Langevin Equations can be combined with a PIMD framework to
accelerate convergence of quantum observables while retaining systematic
approach to the quantum limit. Parameters formatted for i-PI input can
be obtained from the `GLE4MD
website <http://gle4md.org/index.html?page=matrix>`__.

| **Main contributors:** Michele Ceriotti, Joshua More
| **Implementation:**
| M. Ceriotti, J. More, D. Manolopoulos, *“i-PI: A Python interface for
  ab initio path integral molecular dynamics simulations”*, Comp. Phys.
  Comm. 185(3), 1019 (2014)
| DOI:
  `10.1016/j.cpc.2013.10.027 <http://dx.doi.org/10.1016/j.cpc.2013.10.027>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1016%2Fj.cpc.2013.10.027>`__
| **Theory:**
| *PIGLET* — M. Ceriotti and D. E. Manolopoulos, *“Efficient
  First-Principles Calculation of the Quantum Kinetic Energy and
  Momentum Distribution of Nuclei”*, Phys. Rev. Lett. 109, 100604 (2012)
| DOI:
  `10.1103/PhysRevLett.109.100604 <http://dx.doi.org/10.1103/PhysRevLett.109.100604>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1103%2FPhysRevLett.109.100604>`__
| *PI+GLE* — M. Ceriotti, D. E. Manolopoulos, and M. Parrinello,
  *“Accelerating the Convergence of Path Integral Dynamics with a
  Generalized Langevin Equation”*, J. Chem. Phys. 134, 84104 (2011)
| DOI: `10.1063/1.3556661 <http://dx.doi.org/10.1063/1.3556661>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.3556661>`__


Path Integral Molecular Dynamics
--------------------------------

The basic PIMD implementation in i-PI relies on a normal-modes
integrator, and allows setting non-physical masses, so that both RPMD
and CMD can be easily realized.

| **Main contributors:** Michele Ceriotti, Joshua More
| **Implementation:**
| M. Ceriotti, J. More, D. Manolopoulos, *“i-PI: A Python interface for
  ab initio path integral molecular dynamics simulations”*, Comp. Phys.
  Comm. 185(3), 1019 (2014) DOI:
  `10.1016/j.cpc.2013.10.027 <http://dx.doi.org/10.1016/j.cpc.2013.10.027>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1016%2Fj.cpc.2013.10.027>`__
| **Theory:**
| R. Feynman, A. Hibbs, *“Quantum Mechanics and Path Integrals”*,
  McGraw-Hill (1964)
| M. Tuckerman, *“Statistical Mechanics and Molecular Simulations”*,
  Oxford Univ. Press (2008)

Path Integrals Molecular Dynamics at Constant Pressure
------------------------------------------------------

The constant-pressure implementation allows for arbitrary thermostats to
be applied to the cell degrees of freedom, and work in both
constant-shape and variable-cell mode.

| **Main contributors:** Michele Ceriotti, Joshua More, Mariana Rossi
| **Implementation:**
| M. Ceriotti, J. More, D. Manolopoulos, *“i-PI: A Python interface for
  ab initio path integral molecular dynamics simulations”*, Comp. Phys.
  Comm. 185(3), 1019 (2014)
| DOI:
  `10.1016/j.cpc.2013.10.027 <http://dx.doi.org/10.1016/j.cpc.2013.10.027>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1016%2Fj.cpc.2013.10.027>`__
| **Theory:**
| G. J. Martyna, A. Hughes, M. Tuckerman, *“Molecular dynamics
  algorithms for path integrals at constant pressure”*, J. Chem. Phys.
  110(7), 3275 (1999)
| DOI: `10.1063/1.478193 <http://dx.doi.org/10.1063/1.478193>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.478193>`__
| G. Bussi, T. Zykova-Timan, M. Parrinello, *“Isothermal-isobaric
  molecular dynamics using stochastic velocity rescaling”*, J. Chem.
  Phys. 130(7), 074101 (2009)
| DOI: `10.1063/1.3073889 <http://dx.doi.org/10.1063/1.3073889>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.3073889>`__
| P. Raiteri, J. D. Gale, G. Bussi, *“Reactive force field simulation of
  proton diffusion in BaZrO3 using an empirical valence bond approach”*,
  J. Phys. Cond. Matt. 23(33), 334213 (2011)
| DOI:
  `10.1088/0953-8984/23/33/334213 <http://dx.doi.org/10.1088/0953-8984/23/33/334213>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1088/0953-8984/23/33%2F334213>`__

Path Integral Langevin Equation Thermostats
-------------------------------------------

Simple yet efficient Langevin thermostat for PIMD, with normal-modes
thermostats optimally coupled to the ideal ring polymer frequencies

| **Main contributors:** Michele Ceriotti
| **Implementation and Theory:**
| M. Ceriotti, M. Parrinello, T. E. Markland, and D. E. Manolopoulos,
  *“Efficient stochastic thermostatting of path integral molecular
  dynamics”* J. Chem. Phys. 133, 124104 (2010).
| DOI: `10.1063/1.3489925 <http://dx.doi.org/10.1063/1.3489925>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.3489925>`__

Perturbed Path Integrals
------------------------

Effectively a zeroth-order cumulant expansion of the high-order PI
Hamiltonian, perturbed path integrals offer an attractive approach to
compute thermochemistry of materials and molecules including quantum
nuclei, as a post-processing of a Trotter trajectory.

| **Main contributors:** Igor Poltavski
| **Theory:**
| I. Poltavsky and A. Tkatchenko, *“Modeling Quantum Nuclei with
  Perturbed Path Integral Molecular Dynamics”*, Chem. Sci. 7, 1368
  (2016)
| DOI: `10.1039/C5SC03443D <http://dx.doi.org/10.1039/C5SC03443D>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1039%2FC5SC03443D>`__


Replica Exchange MD
-------------------

Accelerated convergence of averages by performing Monte Carlo exchanges
of configurations between parallel calculations.

| **Main contributors:** Riccardo Petraglia, Robert Meissner, Michele
  Ceriotti
| **Implementation:** R. Petraglia, A. Nicolaï, M. M. D. Wodrich, M.
  Ceriotti, C. Corminboeuf, *“Beyond static structures: Putting forth
  remd as a tool to solve problems in computational organic chemistry”*,
  J. Comput. Chem. 37(1), 83-92 (2016)
| DOI: `10.1002/jcc.24025 <http://dx.doi.org/10.1002/jcc.24025>`__ — BibTeX:
  `fetch <http://www.doi2bib.org/#/doi/10.1002/jcc.24025>`__
| **Theory:**
| Y. Sugita, Y. Okamoto, *“Replica-exchange molecular dynamics method
  for protein folding”*, Chem. Phys. Lett. 314(1-2), 141–151 (1999)
| DOI:
  `10.1016/s0009-2614(99)01123-9 <http://dx.doi.org/10.1016/s0009-2614(99)01123-9>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/#/doi/10.1016/s0009-2614(99)01123-9>`__
| T. Okabe, M. Kawata, Y. Okamoto, M. Mikami, *“Replica-exchange monte
  carlo method for the isobaric–isothermal ensemble”*, Chem. Phys. Lett.
  335(5-6), 435-439 (2001)
| DOI: `10.1016/
  s0009-2614(01)00055-0 <http://dx.doi.org/10.1016/s0009-2614(01)00055-0>`__ —
  BibTeX:
  `fetch <http://www.doi2bib.org/#/doi/10.1016/s0009-2614(01)00055-0>`__

Quantum Alchemical Transformation
---------------------------------

An algorithm that performs Monte Carlo moves to change a chemical
species into its isotopes.

| **Main contributors:** Bingqing Cheng, Michele Ceriotti
| **Implementation:**
| Cheng, Bingqing, J"{o}rg Behler, Michele Ceriotti, *“Nuclear Quantum
  Effects in Water at the Triple Point: Using Theory as a Link Between
  Experiments.”* J. Phys. Chem. Lett. 7(12), 2210-2215 (2016)
| DOI:
  `10.1021/acs.jpclett.6b00729 <http://dx.doi.org/10.1021/acs.jpclett.6b00729>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1021%2Facs.jpclett.6b00729>`__
| **Theory:**
| Michael R. Shirts, David L. Mobley, John D. Chodera, *“Alchemical Free
  Energy Calculations: Ready for Prime Time?”*, Ann. Rep. Comp. Chem.
  41-59 (2007)
| DOI:
  `10.1016/S1574-1400(07)03004-6 <http://dx.doi.org/10.1016/S1574-1400(07)03004-6>`__
  — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1016%2FS1574-1400(07)03004-6>`__
| Jian Liu, Richard S Andino, Christina M Miller, Xin Chen, David M
  Wilkins, Michele Ceriotti, David E Manolopoulos, *“A surface-specific
  isotope effect in mixtures of light and heavy water”*, J. Phys. Chem.
  C 117(6), 2944-2951 (2013)
| DOI: `10.1021/jp311986m <http://dx.doi.org/10.1021/jp311986m>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/#/10.1021/jp311986m>`__


Reweighting-based high-order PIMD
---------------------------------

The Boltzmann weight assciated with the high order correction to
standard PIMD is printed out as a property so that the high order
estimate of an arbitrary position-dependent observable can be computed
as a weighted average.

| **Main contributors:** Michele Ceriotti , Guy A. R. Brian
| **Implementation:**
| M.Ceriotti, G.A.R.Brian, O.Riordan, D.E.Manolopolous *“The
  inefficiency of re-weighted sampling and the curse of system size in
  high-order path integration”*, Proc. R. Soc. A 468, 2-17 (2011)
| DOI:
  `10.1098/rspa.2011.0413 <http://dx.doi.org/10.1098/rspa.2011.0413>`__
  — BIBTEX:
  `fetch <http://www.doi2bib.org/bib/10.1098%2Frspa.2011.0413>`__
| **Theory:**
| S.Jang, S.Jang, G.A.Voth *“Applications of higher order composite
  factorization schemes in imaginary time path integral simulations”*,
  J. Chem. Phys. 115, 7832 (2001)
| DOI: `10.1063/1.1410117 <http://dx.doi.org/10.1063/1.1410117>`__ —
  BIBTEX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.1410117>`__
| S.A.Chin *“Symplectic integrators from composite operator
  factorizations”*, Phys. Lett. A 226, 344 (1997)
| DOI:
  `10.1016/s0375-9601(97)00003-0 <http://dx.doi.org/10.1016/s0375-9601(97)00003-0>`__
  — BIBTEX:
  `fetch <http://www.doi2bib.org/bib/10.1016%2Fs0375-9601(97)00003-0>`__
| M.Suzuki *“Hybrid exponential product formulas for unbounded operators
  with possible applications to Monte Carlo simulations”*, Phys. Lett. A
  201, 425 (1995)
| DOI:
  `10.1016/0375-9601(95)00266-6 <http://dx.doi.org/10.1016/0375-9601(95)00266-6>`__
  — BIBTEX:
  `fetch <http://www.doi2bib.org/bib/10.1016%2F0375-9601(95)00266-6>`__

Ring-Polymer Contraction
------------------------

A ring-polymer contraction makes it possible to compute different
components of the forces on different number of imaginary time slices.
In order to reap maximum benefits, the implementation is fully
compatible with the multiple time step integrators.

| **Main contributors:** Michele Ceriotti, Venkat Kapil
| **Implementation:**
| V.Kapil, J.VandeVondele, M.Ceriotti *“Accurate molecular dynamics and
  nuclear quantum effects at low cost by multiple steps in real and
  imaginary time: using density functional theory to accelerate
  wavefunction methods”*, J. Chem. Phys. 144, 054111 (2016) DOI:
  `10.1063/1.4941091 <http://dx.doi.org/10.1063/1.4941091>`__ — BibTeX:
  `fetch <http://www.doi2bib.org/bib/10.1063%2F1.4941091>`__
| **Theory:**
| T.Markland, D.E.Manolopoulos *“An efficient ring polymer contraction
  scheme for imaginary time path integral simulations”*, J. Chem. Phys.
  129, 024105 (2008)
| DOI: `10.1063/1.2953308 <http://dx.doi.org/10.1063/1.2953308>`__ —
  BibTeX: `fetch <http://www.doi2bib.org/bib/10.1063%2F1.2953308>`__


Ring-polymer Instantons
-----------------------

Semiclassical instanton theory is an efficient way of simulating
tunneling contributions to reaction rate constants and tunneling
splittings, based on a well-defined dominant tunneling pathway. It can
be much more efficient than RPMD rate theory, but it is not applicable
to condensed phases and includes anharmonicities only along the reaction
coordinate.

| **Main contributors:** Yair Litman, Jeremy O. Richardson, Mariana
  Rossi
| **Implementation:**
| Y. Litman, J. O. Richardson, T. Kumagai, M. Rossi, *Elucidating the
  Quantum Dynamics of Intramolecular Double Hydrogen Transfer in
  Porphycene*, arXiv:1810.05681 (2018).
| V. Kapil et al. *i-PI 2.0: A Universal Force Engine for
  AdvancedMolecular Simulations*, Comp. Phys. Comm. (2018)
| **Theory:**
| W. H. Miller, *Semiclassical limit of quantum mechanical transition
  state theory for nonseparable systems*, J. Chem. Phys. 62(5) 1899–1906
  (1975)
| DOI: `10.1063/1.430676 <http://dx.doi.org/10.1063/1.430676>`__ — BibTeX:
  `fetch <http://doi2bib.org/#/doi/10.1063/1.430676>`__
| J. O. Richardson, *Ring-polymer instanton theory*, Int. Rev. Phys.
  Chem. 37, 171 (2018)
| DOI:
  `10.1080/0144235X.2018.1472353 <http://doi.org/10.1080/0144235X.2018.1472353>`__
  — BibTeX: `fetch <http://doi2bib.org/#/doi/10.1080/0144235X.2018.1472353>`__

Thermodynamic Integrations
--------------------------

Thermodynamic integrations are made easy with i-PI, through the
connection of different sockets and the different weights one can assign
to different forces. It is possible to do harmonic (Debye model) to
anharmonic integration, or integrations between different kinds of
potentials. Also, quantum thermodynamic integrations relative to mass
are easily done through the manual input of each atom’s masses.

| **Main contributors:** Mariana Rossi, Michele Ceriotti
| **Implementation:**
| M. Rossi, P. Gasparotto, M.Ceriotti *“Anharmonic and Quantum
  Fluctuations in Molecular Crystals: A First-Principles Study of the
  Stability of Paracetamol”*, Phys. Rev. Lett. 117, 115702 (2016)
| DOI:
  `10.1103/PhysRevLett.117.115702 <http://dx.doi.org/10.1103/PhysRevLett.117.115702>`__
  — BIBTEX:
  `fetch <http://www.doi2bib.org/bib/10.1103%2FPhysRevLett.117.115702>`__

Thermostatted Ring-polymer Molecular Dynamics
---------------------------------------------

By introducing an internal mode thermostat to RPMD it is possible to
reduce the well-known artifacts in the simulation of dynamical
properties by path integral methods.

| **Main contributors:** Mariana Rossi, Michele Ceriotti
| **Implementation:**
| M.Rossi, M.Ceriotti, D.E.Manolopoulos, *“How to remove the spurious
  resonances from ring polymer molecular dynamics”*, J. Chem. Phys. 140,
  234116 (2014)
| DOI: `10.1063/1.4883861 <https://doi.org/10.1063/1.4883861>`__ —
  BibTeX: `fetch <https://www.doi2bib.org/bib/10.1063%2F1.4883861>`__
| **Theory:**
| I.R.Craig, D.E.Manolopoulos *“Quantum statistics and classical
  mechanics: Real time correlation functions from ring polymer molecular
  dynamics”*, J. Chem. Phys. 121, 3368 (2004)
| DOI: `10.1063/1.1777575 <https://doi.org/10.1063/1.1777575>`__ —
  BibTeX: `fetch <https://www.doi2bib.org/bib/10.1063%2F1.1777575>`__

