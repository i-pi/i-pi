Open path integral calculations for momentum distribution. 
=================================================

`one_path`: An open-path integral calculation with one atom represented with an open path while all others are represented with a closed path. This setup can be used for rigorous estimation of the momentum distribution of a particle. 

`multi_path`: An open-path integral calculation with multiple atoms of a single atom type are represented with an open path while all others are represented with a closed path. Practically, this allows to average the momentum distribution over multiple atoms but it is not a rigorous way to estimate the momentum distribution of an atom if the kinetic energies of the different atoms are correlated.  

`eco_multi_path`: Same as `multi_path`, but using economised ("eco") path
integrals: the spring frequencies of both the closed and the open paths are
optimized to accelerate convergence with the number of beads, the open ones
being fitted to reproduce the end-to-end distribution (hence the momentum
distribution) and the radius of gyration of open harmonic paths. The
end-to-end vectors are recovered from the difference between the first and
last bead positions (the `posA`/`posB` trajectory outputs).
