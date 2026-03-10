#!/usr/bin/env python3

"""
Computes the autocorrelation function from i-pi outputs.
Assumes the input files are in xyz format and atomic units.
"""


import argparse
import numpy as np

from ipi.utils.messages import verbosity
from ipi.utils.tools import compute_acf_xyz

verbosity.level = "low"

if __name__ == "__main__":
    # adds description of the program.
    parser = argparse.ArgumentParser(
        description="Given a xyz formatted vector, computes its autocorrelation function and its Fourier transform, Parses xyz formatted files with units specified accoridng to i-pi standards. Produces the result in atomic units."
    )

    # adds arguments.
    parser.add_argument(
        "-ifile",
        "--input_file",
        required=True,
        type=str,
        default=None,
        help="the relative path to the xyz formatted file",
    )
    parser.add_argument(
        "-mlag",
        "--maximum_lag",
        required=True,
        type=int,
        default=None,
        help="the maximum time lag for the autocorrelation function",
    )
    parser.add_argument(
        "-bsize",
        "--block_length",
        type=int,
        default=-1,
        help="the number of lines to be imported at once during ``chunk-by-chunk`` input; defaults to 2 * MAXIMUM_LAG",
    )
    parser.add_argument(
        "-ftpad",
        "--length_zeropadding",
        type=int,
        default=0,
        help="number of zeroes to be padded at the end of the autocorrelation function before the Fourier transform is computed",
    )
    parser.add_argument(
        "-ftwin",
        "--spectral_windowing",
        type=str,
        choices=[
            "none",
            "cosine-hanning",
            "cosine-hamming",
            "cosine-blackman",
            "triangle-bartlett",
        ],
        default="none",
        help="type of window function the autocorrelation function is multiplied with before the Fourier transform is computed.",
    )
    parser.add_argument(
        "-dt",
        "--timestep",
        type=str,
        default="1 atomic_unit",
        help="timestep associated with consecutive frames. <number> <unit>. Defaults to 1.0 atomic_unit",
    )
    parser.add_argument(
        "-labels",
        "--labels",
        type=str,
        default="*",
        help="labels of the species to be monitored",
    )
    parser.add_argument(
        "-s",
        "--skip",
        type=int,
        default=0,
        help="number of initial frames to be skipped",
    )
    parser.add_argument(
        "-oprefix",
        "--output_prefix",
        required=True,
        type=str,
        help="the prefix of the output file.",
    )
    parser.add_argument(
        "-der",
        "--derivative",
        action="store_true",
        help="computes the autocorrelation function of the time derivative of the xyz formatted input",
    )

    args = parser.parse_args()
    timestep, time_units = args.timestep.split()
    mlag, npad = args.maximum_lag, args.length_zeropadding

    # calls utility function to evaluate acf and its ft
    time, vvacf, vvacf_err, omega, fvvacf, fvvacf_err = compute_acf_xyz(
        input_file=args.input_file,
        maximum_lag=args.maximum_lag,
        block_length=args.block_length,
        length_zeropadding=args.length_zeropadding,
        spectral_windowing=args.spectral_windowing,
        labels=args.labels.split(),
        timestep=float(timestep),
        time_units=time_units,
        skip=args.skip,
        compute_derivative=args.derivative,
    )

    ofile = args.output_prefix
    np.savetxt(ofile + "_acf.data", np.c_[time, vvacf, vvacf_err][: mlag + npad])
    np.savetxt(ofile + "_facf.data", np.c_[omega, fvvacf, fvvacf_err][: mlag + npad])
