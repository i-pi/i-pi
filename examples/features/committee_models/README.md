Committee models for uncertainty estimation 
===========================================

Committee (ensemble) models provide a simple, general strategy to evaluate the errors in the predictions of a ML model. 
Focusing on the case of the interatomic potential, one evaluates several models in parallel (which differ because of the architecture, the training set, or just the random initialization) using their mean as the best estimate and their standard deviation as the measure of error.

The `<ffcommittee>` forcefield style provides an implementation of the "calibrated committee" framework discussed in [Musil 2019](http://doi.org/10.1021/acs.jctc.8b00959) and [Imbalzano 2021](10.1063/5.0036522). 

All the simulations here use an artificial committee model built by arbitrarily modifying a qTIP4P/f potential. The two resulting forcefields play the role of two members of a committee of ML models. 

`thermodynamic_uq` provides an example of how to run a simulation using `<ffcommittee>`, and to analyze it to compute the propagated uncertainty on the thermodynamic averages. In this example, the various forcefields are run with different drivers.

`weighted_baseline` provides an example of how to use the uncertainty estimate to revert the predictions to a baseline model when the error becomes too large. The ML model in this case is meant to describe a correction over the baseline, which is smoothly suppressed if its uncertainty increases beyond a cutoff.  

`extra_json` shows how to run this kind of simulations with a model that provides all the estimates in one go as part of a JSON formatted `extras`` string

`weighted_baseline_extra` replicates the `weighted_baseline` model using a JSON-based ensemble
