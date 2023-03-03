#!/usr/bin/env python3

import argparse
import numpy as np

description = """
Computes the quantum momentum distribution of a particle given the end-to-end distances.
It computes both the three components of the momentum distribution and the spherically-averaged
distribution of the proton momentum. It also computes <p^2> in the various directions, and the
total contribute.
"""


def r2_K(d, r, is2half):
    """
    Computes the kernel that describes the contribution of an end-to-end distance to the radial average n(r) and multiplies by r squared.
    For small r, uses a truncated Taylor series expansion.

    Parameters
    ----------
    d           :   float
                    The end to end distance.
    r           :   ndarray
                    The grid.
    is2half     :   float
                    Half times inverse sigma squared, where sigma is the width of the 3D Gaussian Kernel.

    Returns
    -------
    d_r2K       :   ndarray
                    r**2 times the radial kernel computed on the entire grid.
    """

    if d <= 1.0e-2:
        r2 = r**2
        d_r2K = (
            np.exp(-is2half * r2) * 4 * r2
            + 4.0
            / 3.0
            * np.exp(-is2half * r2)
            * r2
            * is2half
            * (-3.0 + 2.0 * r2 * is2half)
            * d**2
        )
    else:
        d_r2K = (
            (np.exp(-is2half * (d - r) ** 2) - np.exp(-is2half * (d + r) ** 2))
            * r
            / (d * is2half)
        )
    return d_r2K


def r2_dK(d, r, is2half):
    """
    Computes the kernel that is used for calculating the derivative of the histogram of the end-to-end distance and multiplies by r squared.
    For small r, uses a truncated Taylor series expansion.

    Parameters
    ----------
    d           :   float
                    The end to end distance.
    r           :   ndarray
                    The grid.
    is2half     :   float
                    Half times inverse sigma squared, where sigma is the width of the 3D Gaussian Kernel.

    Returns
    -------
    d_r2_dK     :   ndarray
                    r**2 times the radial kernel computed on the entire grid.
    """

    if d <= 1.0e-2:
        d_r2_dK = (
            (8.0 * r**3 * d)
            * np.exp(-(r**2) * is2half)
            / (3.0 * (1.0 / is2half) ** 2.5 * np.sqrt(np.pi))
        )
    else:
        d_r2_dK = (
            np.sqrt(1.0 / is2half)
            * (
                (1.0 + 2.0 * r * is2half * d) * np.exp(-is2half * (r + d) ** 2)
                + np.exp(-is2half * (r + d) ** 2 + 4 * r * is2half * d)
                * (-1 + 2 * r * is2half * d)
            )
        ) / (2.0 * np.sqrt(np.pi) * d**2)
    return d_r2_dK


def r2_hist_K(qpath, r, r2_kernel, k_param):
    """
    Computes r squared times the histogram of the end-to-end distance.
    It is evaluated  as the sum of r squared times the radial Kernel evaluated on the grid for each instance of the end-to-end distance.

    Parameters
    ----------
    qpath       :   ndarray
                    A numpy array of the position of all the beads of the  open polymer of the target species.
    r           :   ndarray
                    The grid.
    r2_kernel   :   function
                    r squared times the kernel used for binning the end to end distances.
    k_param     :   float
                    Half times inverse sigma squared, where sigma is the width of the corresponding 3D Gaussian Kernel.

    Returns
    -------
    d_r2_hist_K :   ndarray
                    r**2 times the radial kernel computed on the entire grid.
    """

    d_r2_hist_K = np.zeros(r.shape, float) * 0.0
    for q in qpath:
        d = np.linalg.norm(q[:3] - q[-3:])
        d_r2_hist_K += r2_kernel(d, r, k_param)
    return d_r2_hist_K


def dhist_by_dr(qpath, fpath, r, r2_kernel, k_param):
    """
    Computes the derivative of the histogram of the end-to-end distance.
    It is evaluated as the sum of the scaled gradient estimator times a kernel evaluated on the grid for each instance of the imaginary time path q(Tau).

    Parameters
    ----------
    qpath       :   ndarray
                    A numpy array of the position of all the beads of the  open polymer of the target species.
    fpath       :   ndarray
                    A numpy array of the physical force acting on all the beads of the  open polymer of the target species.
    r           :   ndarray
                    The grid.
    r2_kernel   :   function
                    r squared times the kernel used for binning the scaled gradient estimator.
    k_param     :   list
                    The list of parameters of the kernel and the scaled gradient estimator.

    Returns
    -------
    d_dhist_by_dr  :   ndarray
                    r**2 times the radial kernel computed on the entire grid.
    """

    d_dhist_by_dr = np.zeros(r.shape, float) * 0.0

    is2half = k_param[0]
    coeffs = k_param[1]
    beta_P = k_param[2]
    mw_P2 = k_param[3]
    P = qpath.shape[1] / 3

    for i in range(len(qpath)):
        # Instantaneous position of all the beads of the ring polymer
        q = qpath[i]

        # Instantaneous force acting on all the beads of the ring polymer.
        f = fpath[i]

        # Calculates the end-to-end joining vector.
        x = q[:3] - q[-3:]
        d = np.linalg.norm(x)
        spring_q = q.copy()

        # Calculates the derivative of the spring term with the mw_P^2 term.
        spring_q[3 : 3 * (P - 1)] += q[3 : 3 * (P - 1)]
        spring_q[3:] -= q[: 3 * (P - 1)]
        spring_q[: 3 * (P - 1)] -= q[3:]

        # Multiplies the force and the position with the coefficients to calculates the scaled gradient.
        c_f = np.asarray(
            [np.dot(f[0::3], coeffs), np.dot(f[1::3], coeffs), np.dot(f[2::3], coeffs)]
        )
        c_mw_P2_spring_q2 = mw_P2 * np.asarray(
            [
                np.dot(spring_q[0::3], coeffs),
                np.dot(spring_q[1::3], coeffs),
                np.dot(spring_q[2::3], coeffs),
            ]
        )
        scaled_gradient = beta_P * (c_mw_P2_spring_q2 - c_f)

        # Calculates the projection of the scaled gradient on the end-to-end joining vector.
        if d < 1e-3:
            scaled_gradient_dot_x = 0.0
        else:
            scaled_gradient_dot_x = np.dot(scaled_gradient, x) / d

        # Divides the "r**2 times kernel" by r**2 and avoids the 0 / 0 case when r tends to 0.
        d_dhist_by_dr[r > 1e-3] += (scaled_gradient_dot_x * r2_kernel(d, r, is2half))[
            r > 1e-3
        ] / r[r > 1e-3] ** 2

    return d_dhist_by_dr


def get_np(qpath_file, fpath_file, prefix, bsize, P, m, T, s, ns, skip, der):
    """
    Computes the radial distribution of the particle momentum and the end-to-end distance.

    Parameters
    ----------
    qpath_file  :   ndarray
                    A numpy array of the position of all the beads of the open polymer of the target species.
    fpath_file  :   ndarray
                    A numpy array of the physical force acting on all the beads of the open polymer of the target species.
    prefix      :   ndarray
                    The grid.
    block_size  :   function
                    r squared times the kernel used for binning the scaled gradient estimator.
    P           :   float
                    The number of beads of the open polymer.
    m           :   float
                    The mass of the target species (a.m.u.).
    T           :   float
                    The temperature (K).
    s           :   float
                    The extrema of the grid.
    ns          :   int
                    The number of grid points.
    skip        :   int
                    The number of inital `steps' to be skipped.
    der         :   boolean
                    Triggers the scaled derivative estimator.

    """

    # Converts temperature and mass to internal atomic units.
    T = T * 3.1668105e-6
    m = 1822.888 * m
    prefix = prefix + "_"

    # Reads the file containing the position of all the beads of the ring polymer.
    qpath = np.loadtxt(qpath_file, skiprows=int(skip))

    # Reads the file containing the force acting on all the beads of the ring polymer if the flag for calculating the derivative is True.
    if der is True:
        fpath = np.loadtxt(fpath_file, skiprows=int(skip))
        # Defines parameters of the derivative histogram.
        is2half = 0.5 * m * P * T
        coeffs = 0.5 * np.asarray(
            [-1.0 + float(j) * 2.0 / float(P - 1) for j in range(P)]
        )
        beta_P = 1.0 / (P * T)
        mw_P2 = m * (P * T) ** 2
        der_params = [is2half, coeffs, beta_P, mw_P2]

    # Defines the extremum of the grid if not specified.
    if s <= 0:
        s = np.sqrt(np.max(np.sum((qpath[:, :3] - qpath[:, -3:]) ** 2, axis=1))) * 6.0
    # Defines the number of grid points if not specified.
    if ns <= 0:
        ns = int(s * np.sqrt(T * P * m)) * 6 + 1

    # Defines the grids.
    r = np.linspace(0, s, ns)
    dr = abs(r[1] - r[0])
    p = np.linspace(0, np.pi / dr, ns)
    dp = np.abs(p[0] - p[1])

    # Defines empty python lists that store the histograms (or their derivatives) of the end-to-end distance and the momentum.
    r2_4pi_h_list = []
    p2_4pi_np_list = []
    avgp2_list = []
    if der is True:
        dh_by_dr_list = []

    # Calculates the number of chunks for block averaging.
    n_blocks = int(len(qpath) / bsize)

    if n_blocks == 0:
        print("# ERROR: Not enough data to build a block")
        exit()

    for x in range(n_blocks):
        if der is False:
            qpath_block = qpath[x * bsize : (x + 1) * bsize]
            # fpath_block = fpath[x*bsize : (x+1)*bsize]

            # Calculates the radial distribution function of the end-to-end distance.
            r2_4pi_h_block = (
                4.0 * np.pi * r2_hist_K(qpath_block, r, r2_K, 0.5 * T * P * m)
            )
            r2_4pi_h_list.append(r2_4pi_h_block)

            # Calculates the radial distribution of momentum.
            p2_4pi_np_block = np.zeros(p.shape, float)
            # Takes the Sine transform of r2_4pi_h(r) / r with the correct 0/0 limit.
            fsin = np.zeros(ns)
            for i in range(ns):
                fsin[0] = p[i]
                fsin[1:] = np.sin(p[i] * r[1:]) / r[1:]
                p2_4pi_np_block[i] = 4 * np.pi * sum(p[i] * r2_4pi_h_block * fsin) * dr
            p2_4pi_np_list.append(p2_4pi_np_block)

            # Appends the average value of p^2 modulo a normalization
            avgp2_list.append(np.dot(p**2, p2_4pi_np_block) * dp)

        else:
            qpath_block = qpath[x * bsize : (x + 1) * bsize]
            fpath_block = fpath[x * bsize : (x + 1) * bsize]

            # Calculates the derivative of the histogram of the end-to-end distance.
            dh_by_dr_block = dhist_by_dr(qpath_block, fpath_block, r, r2_dK, der_params)
            dh_by_dr_list.append(dh_by_dr_block)

            # Calculates radial distribution of histogram of the end-to-end distance ny integrating the derivative.
            h_block = np.cumsum(
                dhist_by_dr(qpath_block, fpath_block, r, r2_dK, der_params)
            )
            # Applies the boundary condition.
            h_block = h_block - h_block[-1]
            r2_4pi_h_block = 4.0 * np.pi * r**2 * h_block
            r2_4pi_h_list.append(r2_4pi_h_block)

            # Calculates the radial distribution of momentum by an integral by parts.
            p2_4pi_np_block = np.zeros(p.shape, float)
            for i in range(ns):
                r_fcos_pr = r * np.cos(p[i] * r)
                if p[i] < 1e-3:
                    sin_pr_by_p = r
                else:
                    sin_pr_by_p = np.sin(p[i] * r) / p[i]
                # p2_4pi_np_block[i] = 16.0 * np.pi**2 * dr * sum(dh_by_dr_block * (r_fcos_pr - sin_pr_by_p))
                p2_4pi_np_block[i] = dr * sum(
                    dh_by_dr_block * (r_fcos_pr - sin_pr_by_p)
                )
            p2_4pi_np_list.append(p2_4pi_np_block)

            # Appends the average value of p^2 modulo a normalization
            avgp2_list.append(np.dot(p**2, p2_4pi_np_block) * dp)

    # Block averages the radial distribution of the end-to-end distance.
    avg_r2_4pi_h = np.sum(np.asarray(r2_4pi_h_list), axis=0)
    norm_r2_4pi_h = np.sum(avg_r2_4pi_h) * dr
    avg_r2_4pi_h = avg_r2_4pi_h / norm_r2_4pi_h
    err_r2_4pi_h = np.std(
        np.asarray(r2_4pi_h_list) / (norm_r2_4pi_h / n_blocks), axis=0
    ) / np.sqrt(n_blocks)
    np.savetxt(str(prefix + "4pi_r2_h" + ".data"), np.c_[r, avg_r2_4pi_h, err_r2_4pi_h])
    print(
        "# Printing the radial distribution of the end-to-end distance :",
        str(prefix + "4pi_r2_h" + ".data"),
    )

    # Block averages the radial momentum distributions and estimates errors.
    avg_p2_4pi_np = np.sum(np.asarray(p2_4pi_np_list), axis=0)
    norm_p2_4pi_np = np.sum(avg_p2_4pi_np) * dp
    avg_p2_4pi_np = avg_p2_4pi_np / norm_p2_4pi_np
    err_p2_4pi_np = np.std(
        np.asarray(p2_4pi_np_list) / (norm_p2_4pi_np / n_blocks), axis=0
    ) / np.sqrt(n_blocks)
    np.savetxt(
        str(prefix + "4pi_p2_np" + ".data"), np.c_[p, avg_p2_4pi_np, err_p2_4pi_np]
    )
    print(
        "# Printing the radial distribution of the particle momentum :",
        str(prefix + "4pi_p2_np" + ".data"),
    )

    # Also calulates the average value of p^2.
    avg_avgp2 = np.sum(np.asarray(avgp2_list), axis=0) / norm_p2_4pi_np
    err_avgp2 = np.std(
        np.asarray(avgp2_list) / (norm_p2_4pi_np / n_blocks), axis=0
    ) / np.sqrt(n_blocks)
    print("# avg <p2> :", avg_avgp2, "+/-", err_avgp2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-qfile", type=str, help="name of the end-to-end distances file"
    )
    parser.add_argument("-ffile", type=str, default="", help="name of the forces file")
    parser.add_argument(
        "--prefix", type=str, default="out", help="prefix for the output files"
    )
    parser.add_argument(
        "-bsize", type=int, default=50000, help="Specify the size of the blocks"
    )
    parser.add_argument("-P", type=int, default=1, help="Specify the number of beads")
    parser.add_argument(
        "-m",
        type=float,
        default=1.0078,
        help="Specify the mass of the atom in a.m.u. default is hydorgen",
    )
    parser.add_argument(
        "-T",
        type=float,
        default=300,
        help="Specify the temperature of the system in kelvin",
    )
    parser.add_argument(
        "-dint",
        type=float,
        default=0,
        help="Specify the positive extrema of the interval to build the histogram ([-dint,dint])",
    )
    parser.add_argument(
        "-ns",
        type=float,
        default=0,
        help="Specify the number of point to use for the histogram",
    )
    parser.add_argument(
        "-skip", type=int, default=0, help="Specify the number of points to be skipped"
    )
    parser.add_argument(
        "-der",
        action="store_true",
        default=False,
        help="Derives, integrates and then takes the Fourier transform",
    )

    args = parser.parse_args()

    get_np(
        args.qfile,
        args.ffile,
        args.prefix,
        args.bsize,
        args.P,
        args.m,
        args.T,
        args.dint,
        args.ns,
        args.skip,
        args.der,
    )
