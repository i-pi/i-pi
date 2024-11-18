#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from ipi.engine.simulation import Simulation


def direct_reweight(pot, obs, kbT):
    """
    Reweights an observable based on a committee of potentials, at the temperature kT
    ref:
    Torrie, Glenn M., and John P. Valleau. "Nonphysical sampling distributions in Monte
    Carlo free-energy estimation: Umbrella sampling."
    Journal of Computational Physics 23.2 (1977): 187-199.
    https://doi.org/10.1016/0021-9991(77)90121-8

    Inputs:
    pot          : an array (nframes x ncommittee) containing the values of the potentials for each committee
    obs          : an array containing the values of the observables
    kbT          : the temperature, in energy units

    Returns:
    obs_avg_rew  : the observable reweighted for each potential
    weights      : the weights computed for each model and frame
    """

    beta = 1.0 / kbT
    num_pot_frames = pot.shape[0]
    num_obs_frames = obs.shape[0]
    if num_pot_frames != num_obs_frames:
        raise RuntimeError(
            "Potential and observable files have different numbers of frames"
        )

    num_pot_models = pot.shape[1]
    if obs.ndim == 1:
        obs = np.expand_dims(obs, axis=1)

    # vectorized evaluation of weights
    weights = np.array(-beta * (pot.T - np.mean(pot, axis=1)).T, dtype=np.float128)
    weights = np.exp(weights)

    # 1-d array of lenght num_pot_models: normalization of the weights over frames
    norm = np.sum(weights, axis=0)
    for i in range(num_pot_models):
        weights[:, i] /= norm[i]
    obs_avg_rew = obs.T @ weights

    return obs_avg_rew, weights


def CEA(pot, obs, kbT):
    """
    Reweights the quantity by cumulant expansion approximation (CEA)
    ref:
    Michele Ceriotti,  Guy A. R. Brain, Oliver Riordan and David E.
    Manolopoulos 2011 The inefficiency of re-weighted sampling and
    the curse of system size in high-order path integration
    Proc. R. Soc. A.4682–17
    http://doi.org/10.1098/rspa.2011.0413

    Inputs:
    pot          : an array (nframes x ncommittee) containing the values of the potentials for each committee
    obs          : an array containing the values of the observables
    kbT          : the temperature, in energy units

    Returns:
    obs_avg_rew  : the observable reweighted for each potential
    weights      : the weights computed for each model and frame

    """
    beta = 1.0 / kbT
    num_pot_frames = pot.shape[0]
    num_obs_frames = obs.shape[0]
    if num_pot_frames != num_obs_frames:
        raise RuntimeError(
            "Potential and observable files have different numbers of frames"
        )
    num_pot_models = pot.shape[1]
    if obs.ndim == 1:
        obs = np.expand_dims(obs, axis=1)
    num_obs_models = obs.shape[1]

    # get the h matrix h = beta(V^i - barV)
    h_matrix = np.zeros((num_pot_frames, num_pot_models))
    # fast way avoiding double loop
    h_matrix = beta * (pot.T - np.mean(pot, axis=1)).T
    # get frame-average of h, h_avg
    h_avg = np.mean(h_matrix, axis=0)  # 1d-array of length num_pot_models
    # get frame-average of observable, obs_avg
    obs_avg = np.mean(obs, axis=0)  # 1d-array of length num_obs_models
    # get frame-average of the product between the observable and h, obsxh_avg
    obsxh_avg = (obs.T @ h_matrix) / num_pot_frames
    # get frame-average of observable, reweighted according to the CEA, obs_avg_CEA
    obs_avg_CEA = np.zeros((num_obs_models, num_pot_models))
    for j in range(num_obs_models):
        for i in range(num_pot_models):
            obs_avg_CEA[j, i] = obs_avg[j] - obsxh_avg[j, i] + obs_avg[j] * h_avg[i]
    return obs_avg_CEA, h_matrix


def uncertainty_CEA_multiple_models(pot, obs, kbT):
    """
    Reweights the quantity by cumulant expansion approximation (CEA) in the case of
    multiple models used to predict the observable.
    ref:
    Michele Ceriotti,  Guy A. R. Brain, Oliver Riordan and David E.
    Manolopoulos 2011 The inefficiency of re-weighted sampling and
    the curse of system size in high-order path integration
    Proc. R. Soc. A.4682–17
    http://doi.org/10.1098/rspa.2011.0413

    Inputs:
    pot         : an array (nframes x ncommittee) containing the values of the potentials for each committee
    obs         : an array containing the values of the observables
    kbT         : the temperature, in energy units

    Returns:
    mean_value  : the average value of the observable
    sigma2_a    : the uncertainty related to the models used to predict the observable
    sigma2_aV   : the uncertainty related to the use of a committee of potentials driving the dynamics
    sigma2_tile : the total uncertainty for the observable

    """
    beta = 1.0 / kbT  # noqa
    obs_avg = np.mean(obs, axis=0)
    obs_avg_CEA, _h_matrix = CEA(obs, pot, kbT)
    fac_a = (pot.shape[1] * (obs.shape[1] - 1)) / (pot.shape[1] * obs.shape[1] - 1)
    fac_aV = ((pot.shape[1] - 1) * obs.shape[1]) / (pot.shape[1] * obs.shape[1] - 1)
    sigma2_a = np.var(obs_avg, ddof=1)  # Eq.(27)
    sigma2_aV = np.mean(np.var(obs_avg_CEA, axis=1, ddof=1))  # Eq.(28)
    sigma2_tilde = fac_a * sigma2_a + fac_aV * sigma2_aV  # Eq.(26)
    return np.mean(obs_avg), sigma2_a, sigma2_aV, sigma2_tilde


def commitee_reweight(
    pot_file, obs_file, kt, stride=1, index=-1, direct=False, multi_models=False
):
    """
    Parameters
    ----------
    pot_file    :   str, mandatory
                    The file containing the values of the potential from each model, either in .dat form or as .extra.
                    Each line contains the potentials for the M models for each frame
    obs_file    :   str, mandatory
                    A file containing the value of the observable for each frame.
                    It is assumed that lines correspond to frames, while columns to properties.
                    Multiple properties can be reweighted at the same time.
    kt          :   float, mandatory
                    The thermal energy, in the same units as pot_file
    stride      :   integer, [1]
                    The frequency of sampling of pot_file, if different from that of prop_file (e.g. --stride 10 if
                    the potential was printed 10 times more often than the observable. --stride -10 must be used if
                    the observable were printed more often than the potential).
    index       :   integer, optional
                    Which column of the property file should be used to compute the average and uncertainty.
                    If it is not given, and there are multiple columns, they will be interpreted as corresponding
                    to multiple property models
    direct      :   bool, optional
                    Activates the direct reweighting, instead of the cumulant expansion approximation.
                    Use at your own risk! Does not work with multi_models
    multi_models:   bool, optional
                    Activates the uncertainty for multiple models, as in Eqs. 26,27,28 of the paper.
                    It considers each column as a prediction from a single model and returns the
                    uncertainty of the models together with the uncertainty of the potentials.

    """
    if index >= 0:
        obs = np.loadtxt(obs_file, usecols=index)
    else:
        obs = np.loadtxt(obs_file)

    if stride > 0:
        potentials = np.loadtxt(pot_file)[::stride]
    elif stride < 0:
        stride = np.abs(stride)
        obs = obs[::stride]
    else:
        raise ValueError("Stride value cannot be zero")

    if multi_models:
        mean_value, sigma2_a, sigma2_aV, sigma2_tilde = uncertainty_CEA_multiple_models(
            potentials, obs, kt
        )
        print("#     Mean            Error         sigma_a        sigma_aV")
        print(
            "{:.8e}  {:.8e}  {:.8e}  {:.8e}".format(
                mean_value, np.sqrt(sigma2_tilde), np.sqrt(sigma2_a), np.sqrt(sigma2_aV)
            )
        )
    else:
        # CEA is the default choice. The weights or h_matrix are
        if direct:
            rw_obs, _weights = direct_reweight(potentials, obs, kt)
        else:
            rw_obs, _h_matrix = CEA(potentials, obs, kt)
        print(
            "#     Mean            Error        <committee_1>         ....         <committee_N>"
        )

        np.savetxt(
            sys.stdout,
            np.vstack([rw_obs.mean(axis=1), rw_obs.std(axis=1), rw_obs.T]).T,
            fmt="% 15.8e ",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This tool exploits a committee of ML models to compute the reweighted  canonical averages of a given observable, in order to quantify its uncertainty at thermodynamic equilibrium.
    The units for the potentials are assumed to be Hartree.
    The full methodology is described in: https://doi.org/10.1063/5.0036522"""
    )

    parser.add_argument(
        "pot_file",
        type=str,
        help="The file containing the potentials, in units of kT (default: Hartree). Rows = frames, columns = potentials",
    )
    parser.add_argument(
        "obs_file",
        type=str,
        help="The file containing the properties. Rows = frames, columns = property/ies",
    )
    parser.add_argument(
        "--kt",
        type=float,
        default=0.0,
        help="The thermal energy. Should be in the same units of the potentials, typically Hartree",
    )
    parser.add_argument(
        "--input",
        type=str,
        default="",
        help="The path to the input file used to run the simulation (usually input.xml). Used just to extract kT.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="The stride ratio used to print the potential and the property. A positive number will take one every n-th value of the potential. A negative number will do the same on the property file.",
    )
    parser.add_argument(
        "--index",
        type=int,
        default=-1,
        help="0-based index of the property that we want to compute",
    )
    parser.add_argument(
        "--direct",
        action="store_true",
        help="Call this option to activate direct reweighting. Standard reweighting is the CEA",
    )
    parser.add_argument(
        "--multi-obs",
        action="store_true",
        help="Call this option to activate reweighting for multiple observable models. The obs_file is interpreted as having one row per step, and columns correspond to the members of a property committee. It returns both the uncertainty associated with the spread of the observable models, and the one derived from the reweighting due to the potential models",
    )

    args = parser.parse_args()

    # Load kbT from i-PI, we could make it into a small function
    if args.input != "":
        simul = Simulation.load_from_xml(
            open(args.input), custom_verbosity="quiet", read_only=True
        )
        kt = float(simul.syslist[0].ensemble.temp)
    else:
        kt = args.kt
    if kt <= 0:
        raise ValueError("Must specify either --kt or --input")

    sys.exit(
        commitee_reweight(
            args.pot_file,
            args.obs_file,
            kt,
            args.stride,
            args.index,
            args.direct,
            args.multi_obs,
        )
    )
