import numpy as np
from ipi.utils.io import read_file_raw
from ipi.utils.units import unit_to_internal


def compute_acf_xyz(
    input_file,
    maximum_lag,
    block_length=None,
    length_zeropadding=0,
    spectral_windowing="none",
    labels=["*"],
    atom_mask=None,
    timestep=1.0,
    time_units="atomic_unit",
    skip=0,
    compute_derivative=False,
):
    """
    Helper functions to compute autocorrelation functions (and their transform)
    from an xyz-formatted input (so these are "vectorial" ACFs summed over all
    atoms and Cartesian coordinates).

    :param input_file: name of the file holding the xyz-formatted trajectory
    :param maximum_lag: maximum lag time to compute ACF for
    :param block_length: length of trajectory blocks used for averaging
    :param length_zeropadding: pad with these many zeros to increase resolution
    :param spectral_windowing: window function to smoothen FT
    :param labels: list of (space separated) atom labels to include, "*" for all
    :param atom_mask: None, or a list of 0 and 1 to act as masks for selecting atoms
    :param timestep: the time step between two frames
    :param time_units: string providing the units of time given
    :param skip: how many frames to skip before starting accumulating ACF
    :param compute_derivative: bool (default false) whether to compute derivative ACF
    """
    # stores the arguments
    ifile = input_file
    mlag = maximum_lag or 100
    bsize = block_length or (2 * mlag + 1)
    npad = length_zeropadding
    ftbox = spectral_windowing
    fskip = skip
    der = compute_derivative

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

    # stores the indices of the "chosen" atoms.
    ndof = len(rr["data"])
    if atom_mask is None:
        if "*" in labels:
            labelbool = np.ones(ndof // 3, bool)
        else:
            labelbool = np.zeros(ndof // 3, bool)
            for l in labels:
                labelbool = np.logical_or(labelbool, rr["names"] == l)
    else:
        labelbool = atom_mask
    nspecies = labelbool.sum()

    # initializes variables.
    nblocks = 0
    dt = unit_to_internal("time", time_units, timestep)
    data = np.zeros((bsize, nspecies, 3), float)
    time = np.asarray(list(range(mlag + 1))) * dt
    omega = (
        np.asarray(list(range(2 * (mlag + npad))))
        / float(2 * (mlag + npad))
        * (2 * np.pi / dt)
    )
    fvvacf = np.zeros_like(omega)
    fvvacf2 = np.zeros_like(omega)
    vvacf = np.zeros_like(time)
    vvacf2 = np.zeros_like(time)

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
            fvvacf2 += mfpvvacf**2
            vvacf += mvvacf
            vvacf2 += mvvacf**2

            nblocks += 1

        except EOFError:
            break
    ff.close()

    # Performs the block average of the Fourier transform.
    fvvacf = fvvacf / nblocks
    fvvacf_err = np.sqrt(fvvacf2 / nblocks - fvvacf**2)

    # Computes the inverse Fourier transform to get the vvac.
    vvacf = vvacf / nblocks
    vvacf_err = np.sqrt(vvacf2 / nblocks - vvacf**2)

    return time, vvacf, vvacf_err, omega, fvvacf, fvvacf_err
