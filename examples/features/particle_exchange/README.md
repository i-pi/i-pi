Exchange in path integral molecular dynamics. 
=================================================
Path integral molecular dynamics including bosonic exchange effects (ordinary path integral molecular dynamics assumes that particles are distinguishable). Exchange manifests in ring polymers that link several particles together.
Fermionic statistics can also be obtained by reweighting bosonic simulations, when the sign problem is not too large (see Hirshberg et al., doi:10.1063/5.0008720).
The algorithm is based on the evaluation of an effective spring potential that averages over multiple path connectivities, that scales quadratically with the number of particles and linearly with the number of beads, and is based on Hirshberg et al.'s doi:10.1073/pnas.1913365116 and Feldman and Hirshberg's doi:10.1063/5.0173749.

`basic_trapped_bosons`: trapped bosons and a mixture of bosons and distinguishable particles
`specify_by_label_and_estimators`: example of specifying bosons and kinetic energy estimators via labels
`fermions`: obtaining fermionic statistics by reweighting bosonic simulations
`exchange_probabilities`: properties that indicate the extent that exchange interaction is occurring in the simulation
