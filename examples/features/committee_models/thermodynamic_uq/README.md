Thermodynamic uncertainty estimation
====================================

This runs a simulation in which several potentials are computed, and 
the uncertainty in the potential correction is estimated from the spread of the
members of the comittee, as in [Imbalzano 2021](http://doi.org/10.1021/acs.jctc.8b00959).

Here we take as a simple example that does not require specialized ML software
a simulation of water in which the members differ by the application of a small 
LJ-like perturbation.
The simulation prints out files that are needed to compute how the error in the 
potential propagates on thermodynamic properties, specifically on the g(r).

This also showcases the "active learning" feature, in which structures
with a large error are printed out so they can be further analyzed.

The reweight tool is provided as tools/py/reweight.py and it can used from 
command line. Additionally, the functions can be imported to be run in a 
interactive manner. After running the simulation hereby provided, one can do 
the reweighting and thus access the uncertainty by running the tool as

```
i-pi-committee-reweight input.xml simulation.committee_pot_0 simulation.out --stride 2
```

The `--stride` flag is necessary because the `simulation.committee_pot_0`
has been printed every 5 steps, while the `simulation.out` file has 
been printed every 10 steps.

This will create a file named `rw_simulation.out` where every line is the
reweighted observable computed for every member of the committee. In this
case all of the properties have been reweighted (although it does not make 
sense to reweight the steps).

If one wants to narrow the reweighting to a single property, then they can
add the `--index` flag to indicate the property that must be reweighted.
To reweight the volume, for example, one should use

```
i-pi-committee-reweight input.xml simulation.committee_pot_0 simulation.out --stride 2 --index 5
```

To compute uncertainties in the case in which one also obtains multiple 
predictions for the observables, one can run with the `--multi-obs` option.
For instance, to compute the error on the potential energy (which itself 
comes with multiple predictions) one can run 

```
i-pi-committee-reweight input.xml simulation.committee_pot_0 simulation.committee_pot_0 --multi-obs
```

