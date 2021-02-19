#!/usr/bin/python
import sys
import argparse
import numpy as np
from ipi.inputs.simulation import InputSimulation
from ipi.utils.io.inputs import io_xml


def direct_reweight(pot, obs, kbT):
    """
    reweights an observable based on a committee of potentials, at the temperature kT

    pot: an array (nframes x ncommittee) containing the values of the potentials for each committee
    obs: an array containing the values of the observables
    kbT: the temperature, in energy units
    """
    beta = 1.0 / kbT
    num_pot_frames = pot.shape[0]
    num_obs_frames = obs.shape[0]
    if num_pot_frames != num_obs_frames:
        raise RuntimeError(
            "potential and observable files have different numbers of frames"
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
    reweight the quantity by cumulant expansion approximation (CEA)
    ref:
    Michele Ceriotti,  Guy A. R. Brain, Oliver Riordan and David E.
    Manolopoulos 2011The inefficiency of re-weighted sampling and
    the curse of system size in high-order path integration
    Proc. R. Soc. A.4682–17
    http://doi.org/10.1098/rspa.2011.0413
    :return:
    """
    beta = 1.0 / kbT
    num_pot_frames = pot.shape[0]
    num_obs_frames = obs.shape[0]
    if num_pot_frames != num_obs_frames:
        print("potential and observable have different number of frames")
        sys.exit("exiting...")
    num_pot_models = pot.shape[1]
    if obs.ndim == 1:
        obs = np.expand_dims(obs, axis=1)
    num_obs_models = obs.shape[1]

    print("Computing h matrix")
    # get the h matrix h = beta(V^i - barV)
    h_matrix = np.zeros((num_pot_frames, num_pot_models))
    # fast way avoiding double loop
    h_matrix = beta * (pot.T - np.mean(pot, axis=1)).T
    print("Computing averages")
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


def commitee_reweight(path2ixml, pot_file, obs_file, stride=1, index=-1, direct=False):
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
    stride      :   integer [1]
                    The frequency of sampling of prop_file, if different from that of pot_file (e.g. --stride 10 if
                    observables are output 10 times more rarely than potential
    index       :   integer, optional
                    Which column of the property file should be used to compute the average and uncertainty.
                    If it is not given, and there are multiple columns, they will be interpreted as corresponding
                    to multiple property models
    direct      :   bool, optional
                    Activates the direct reweighting, instead of the cumulant expansion approximation.
                    Use at your own risk!

    """
    potentials = np.loadtxt(pot_file)

    if index >= 0:
        obs = np.loadtxt(obs_file, usecols=index)[::stride]
    else:
        obs = np.loadtxt(obs_file)[::stride]

    # Load kbT from i-PI, we could make it into a small function
    ifile = open(path2ixml, "r")
    xmlrestart = io_xml.xml_parse_file(ifile)
    ifile.close()

    isimul = InputSimulation()
    isimul.parse(xmlrestart.fields[0][1])
    simul = isimul.fetch()

    kbT = float(simul.syslist[0].ensemble.temp)
    # CEA is the default choice. The weights or h_matrix are
    if direct:
        rw_obs, _weights = direct_reweight(potentials, obs, kbT)
    else:
        rw_obs, _h_matrix = CEA(potentials, obs, kbT)

    np.savetxt("rw_" + obs_file + ".dat", rw_obs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""This tool exploits a committee of ML models to compute the reweighted  canonical averages of a given observable, in order to quantify its uncertainty at thermodynamic equilibrium.
    The units are assumed to be atomic, which are also the internal units of i-PI.
    The full methodology is described in: https://doi.org/10.1063/5.0036522"""
    )

    parser.add_argument(
        "input_xml",
        type=str,
        help="The path to the input file used to run the simulation (usually input.xml)",
    )
    parser.add_argument(
        "pot_file",
        type=str,
        help="The file containing the potentials. Rows = frames, columns = potentials",
    )
    parser.add_argument(
        "obs_file",
        type=str,
        help="The file containing the properties. Rows = frames, columns = property/ies",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="The difference in stride between the potential and properties",
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

    args = parser.parse_args()
    sys.exit(
        commitee_reweight(
            args.input_xml,
            args.pot_file,
            args.obs_file,
            args.stride,
            args.index,
            args.direct,
        )
    )
