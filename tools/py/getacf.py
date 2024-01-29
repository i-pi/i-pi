#!/usr/bin/env python3

"""
Computes the autocorrelation function from i-pi outputs. Assumes the input files are in xyz format and atomic units.
"""


import argparse
import numpy as np
from ipi.utils.io import read_file_raw
from ipi.utils.units import unit_to_internal
from ipi.utils.messages import verbosity

verbosity.level = "low"


def compute_acf(
    input_file,
    output_prefix,
    maximum_lag,
    block_length,
    length_zeropadding,
    spectral_windowing,
    labels,
    timestep,
    skip,
    der,
):
    # stores the arguments
    ifile = str(input_file)
    ofile = str(output_prefix)
    mlag = int(maximum_lag)
    bsize = int(block_length)
    npad = int(length_zeropadding)
    ftbox = str(spectral_windowing)
    labels = str(labels).split()
    timestep = str(timestep).split()
    fskip = int(skip)

    # checks for errors
    if mlag <= 0:
        raise ValueError("MAXIMUM_LAG should be a non-negative integer.")
    if npad < 0:
        raise ValueError("LENGTH_ZEROPADDING should be a non-negative integer.")
    if bsize < 2 * mlag:
        if bsize == -1:
            bsize = 2 * mlag
        else:
            raise ValueError(
                "LENGTH_BLOCK should be greater than or equal to 2 * MAXIMUM_LAG."
            )

    # reads one frame.
    ff = open(ifile)
    rr = read_file_raw("xyz", ff)
    ff.close()

    # appends "der" to output file in case the acf of the derivative is desired
    if der is True:
        ofile = ofile + "_der"

    # stores the indices of the "chosen" atoms.
    ndof = len(rr["data"])
    if "*" in labels:
        labelbool = np.ones(ndof // 3, bool)
    else:
        labelbool = np.zeros(ndof // 3, bool)
        for l in labels:
            labelbool = np.logical_or(labelbool, rr["names"] == l)
    nspecies = labelbool.sum()

    # initializes variables.
    nblocks = 0
    dt = unit_to_internal("time", timestep[1], float(timestep[0]))
    data = np.zeros((bsize, nspecies, 3), float)
    time = np.asarray(list(range(mlag + 1))) * dt
    omega = (
        np.asarray(list(range(2 * (mlag + npad))))
        / float(2 * (mlag + npad))
        * (2 * np.pi / dt)
    )
    fvvacf = omega.copy() * 0.0
    fvvacf2 = fvvacf.copy() * 0.0
    vvacf = time.copy() * 0.0
    vvacf2 = time.copy() * 0.0

    # selects window function for fft.
    if ftbox == "none":
        win = np.ones(2 * mlag + 1, float)
    elif ftbox == "cosine-hanning":
        win = np.hanning(2 * mlag + 1)
    elif ftbox == "cosine-hamming":
        win = np.hamming(2 * mlag + 1)
    elif ftbox == "cosine-blackman":
        win = np.blackman(2 * mlag + 1)
    elif ftbox == "triangle-bartlett":
        win = np.bartlett(2 * mlag + 1)

    ff = open(ifile)
    # Skips the first fskip frames
    for x in range(fskip):
        rr = read_file_raw("xyz", ff)

    while True:
        try:
            # Reads the data in blocks.
            for i in range(bsize):
                rr = read_file_raw("xyz", ff)
                data[i] = rr["data"].reshape((ndof // 3, 3))[labelbool]

            if der is True:
                data = np.gradient(data, axis=0) / dt

            # Computes the Fourier transform of the data.
            fdata = np.fft.rfft(data, axis=0)

            # Computes the Fourier transform of the vvac applying the convolution theorem.
            tfvvacf = fdata * np.conjugate(fdata)

            # Averages over all species and sums over the x,y,z directions. Also multiplies with the time step and a prefactor of (2pi)^-1.
            mfvvacf = (
                3.0 * np.real(np.mean(tfvvacf, axis=(1, 2))) * dt / (2 * np.pi) / bsize
            )

            # Computes the inverse Fourier transform to get the vvac.
            mvvacf = np.fft.irfft(mfvvacf)[: mlag + 1]

            # Applies window in one direction and pads the vvac with zeroes.
            mpvvacf = np.append(mvvacf * win[mlag:], np.zeros(npad))

            # Recomputes the Fourier transform assuming the data is an even function of time.
            mfpvvacf = np.fft.hfft(mpvvacf)

            # Accumulates the (f)acfs and their squares.
            fvvacf += mfpvvacf
            fvvacf2 += mfpvvacf ** 2
            vvacf += mvvacf
            vvacf2 += mvvacf ** 2

            nblocks += 1

        except EOFError:
            break
    ff.close()

    # Performs the block average of the Fourier transform.
    fvvacf = fvvacf / nblocks
    fvvacf_err = np.sqrt(fvvacf2 / nblocks - fvvacf ** 2)

    np.savetxt(ofile + "_facf.data", np.c_[omega, fvvacf, fvvacf_err][: mlag + npad])

    # Computes the inverse Fourier transform to get the vvac.
    vvacf = vvacf / nblocks
    vvacf_err = np.sqrt(vvacf2 / nblocks - vvacf ** 2)
    np.savetxt(ofile + "_acf.data", np.c_[time, vvacf, vvacf_err][: mlag + npad])


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

    # Process everything.
    compute_acf(
        args.input_file,
        args.output_prefix,
        args.maximum_lag,
        args.block_length,
        args.length_zeropadding,
        args.spectral_windowing,
        args.labels,
        args.timestep,
        args.skip,
        args.derivative,
    )
