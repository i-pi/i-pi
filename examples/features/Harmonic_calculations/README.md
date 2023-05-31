Harmonic Calculations
=================================================

FD-Phonons: Displaces every degree of freedom by +/- epsilon. Performs finite-difference of forces to estimate the (3N-3, 3N-3) Hessian and estimates phonons for a solid with N atoms.
FD-NormalModes: Displaces every degree of freedom by +/- epsilon.  Performs finite-difference of forces to estimate the (3N-5/6, 3N-5/6) Hessian and estimates the normalmodes of a molecule with N atoms.
eNMFD-NormalModes: Requires a trail Hessian as input. Displaces every normal mode coordinate by +/- epsilon, where epsilon is calculated from an input energy estimate. Performs finite-difference of forces to estimate the (3N-5/6, 3N-5/6) Hessian and estimates the normalmodes of a molecule with N atoms.
