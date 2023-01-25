#!/usr/bin/env python3
"""
Computes the quantum momentum distribution of a particle given the end-to-end distances.
It computes both the three components of the momentum distribution and the radial function.
Moreover it computes <p^2> both in the various directions and the total contribute.
"""

import argparse
import numpy as np
import time


def kernel(x, invsigma=1.0):
    return np.exp(-0.5 * (x * invsigma) ** 2)


def histo_der(qdata, fdata, grid, k, invsigma, m, P, T):
    ly = grid * 0.0
    ns = len(ly)
    dx = grid[1] - grid[0]
    dj = int(8 * invsigma / dx)

    c = 0.5 * np.asarray([-1.0 + float(j) * 2.0 / float(P - 1) for j in range(P)])
    bp = 1.0 / (P * T)
    mwp2 = m * (P * T) ** 2
    for i in range(len(qdata)):
        q = qdata[i]
        f = fdata[i]
        x = q[0] - q[-1]
        jx = int(x / dx + ns / 2.0)
        fc = (f * c).sum()
        s = q.copy()
        s[1 : P - 1] += q[1 : P - 1]
        s[1:] -= q[: P - 1]
        s[: P - 1] -= q[1:]
        sc = -mwp2 * (s * c).sum()
        ly[jx - dj : jx + dj + 1] += (
            -bp * (fc + sc) * k(grid[jx - dj : jx + dj + 1] - x, invsigma)
        )
    return ly * np.sqrt(0.5 / np.pi * invsigma**2)


def histo(data, grid, k, invsigma):
    ly = grid * 0.0
    ns = len(ly)
    dx = grid[1] - grid[0]
    dj = int(8 * invsigma / dx)  # number of standard deviations to be computed
    for i in range(len(data)):
        x = data[i]
        jx = int(x / dx + ns / 2.0)
        ly[jx - dj : jx + dj + 1] += k(grid[jx - dj : jx + dj + 1] - x, invsigma)
    return ly * np.sqrt(1.0 / 2.0 / np.pi * invsigma**2)


def get_np(qfile, ffile, prefix, bsize, P, mamu, Tkelv, s, ns, der, skip):
    start = time.clock()
    prefix = prefix + "_"

    T = Tkelv * 3.1668105e-6
    m = 1822.888 * mamu

    # initialises grids.
    dq = np.zeros((bsize, 3), float)
    ns = int(float(ns) / 2) * 2 + 1
    dqxgrid = np.linspace(-s, s, ns)
    dqygrid = np.linspace(-s, s, ns)
    dqzgrid = np.linspace(-s, s, ns)
    dqxstep = dqxgrid[1] - dqxgrid[0]
    dqystep = dqygrid[1] - dqygrid[0]
    dqzstep = dqzgrid[1] - dqzgrid[0]

    hxlist = []
    hylist = []
    hzlist = []

    nplistx = []
    nplisty = []
    nplistz = []

    # Defines the grid for momentum.
    pxi = -np.pi / (dqxgrid[1] - dqxgrid[0])
    pxf = +np.pi / (dqxgrid[1] - dqxgrid[0])
    pxgrid = np.linspace(pxi, pxf, ns)
    pxstep = np.abs(pxgrid[1] - pxgrid[0])

    pyi = -np.pi / (dqygrid[1] - dqygrid[0])
    pyf = +np.pi / (dqygrid[1] - dqygrid[0])
    pygrid = np.linspace(pyi, pyf, ns)
    pystep = np.abs(pygrid[1] - pygrid[0])

    pzi = -np.pi / (dqzgrid[1] - dqzgrid[0])
    pzf = +np.pi / (dqzgrid[1] - dqzgrid[0])
    pzgrid = np.linspace(pzi, pzf, ns)
    pzstep = np.abs(pzgrid[1] - pzgrid[0])

    if ns % 2 == 0:
        pxgrid = pxgrid - pxstep / 2.0
        pygrid = pygrid - pystep / 2.0
        pzgrid = pzgrid - pzstep / 2.0

    # Read the end to end distances from file
    bq = np.loadtxt(qfile)[int(skip) :]
    if der is True:
        bf = np.loadtxt(ffile)[int(skip) :]

    step = np.shape(bq)[0]

    n_block = int(step / bsize)

    if n_block == 0:
        print("not enough data to build a block")
        exit()

    if der is False:
        for x in range(n_block):
            print("# building the histogram for block $", x + 1)
            dq = bq[x * bsize : (x + 1) * bsize]
            hx = histo(
                dq[:, 0] - dq[:, 3 * (P - 1)], dqxgrid, kernel, np.sqrt(T * P * m)
            )
            hx = (hx + hx[::-1]) * 0.5
            hy = histo(
                dq[:, 1] - dq[:, 3 * (P - 1) + 1], dqygrid, kernel, np.sqrt(T * P * m)
            )
            hy = (hy + hy[::-1]) * 0.5
            hz = histo(
                dq[:, 2] - dq[:, 3 * (P - 1) + 2], dqzgrid, kernel, np.sqrt(T * P * m)
            )
            hz = (hz + hz[::-1]) * 0.5
            hxlist.append(hx)
            hylist.append(hy)
            hzlist.append(hz)

            # Computes the Fourier transform of the end to end vector.
            npx = hx * 0.0
            npy = hy * 0.0
            npz = hz * 0.0
            print("# computing FT for block #", x + 1)
            for i in range(len(hx)):
                npx[i] = (hx * np.cos(pxgrid[i] * dqxgrid)).sum() * dqxstep
                npy[i] = (hy * np.cos(pygrid[i] * dqygrid)).sum() * dqystep
                npz[i] = (hz * np.cos(pzgrid[i] * dqzgrid)).sum() * dqzstep

            nplistx.append(npx)
            nplisty.append(npy)
            nplistz.append(npz)
    else:
        for x in range(n_block):
            print("# building the histogram for block $", x + 1)
            dq = bq[x * bsize : (x + 1) * bsize]
            df = bf[x * bsize : (x + 1) * bsize]

            c = 0.5 * np.asarray([-1 + 2 * float(j) / float(P - 1) for j in range(P)])
            mwp2 = m * (P * T) ** 2 / (P - 1)
            bp = 1.0 / (P * T)
            for i in range(len(dq)):
                q = dq[i]
                f = df[i]
                x = q[:3] - q[-3:]
                fc = np.asarray(
                    [(f[0::3] * c).sum(), (f[1::3] * c).sum(), (f[2::3] * c).sum()]
                )
                s = q.copy()
                s[3 : 3 * (P - 1)] += q[3 : 3 * (P - 1)]
                s[3:] -= q[: 3 * (P - 1)]
                s[: 3 * (P - 1)] -= q[3:]
                sc = np.asarray(
                    [(s[0::3] * c).sum(), (s[1::3] * c).sum(), (s[2::3] * c).sum()]
                )
                sc *= -mwp2
                fi = -bp * (fc + sc)
                print("@*@", x, fi)
                print("@@@", np.sqrt((x * x).sum()), np.dot(x, fi))

            hx = histo_der(
                dq[:, 0::3], df[:, 0::3], dqxgrid, kernel, np.sqrt(T * P * m), m, P, T
            )
            hx = np.cumsum((hx - hx[::-1]) * 0.5) * dqxstep
            hy = histo_der(
                dq[:, 1::3], df[:, 1::3], dqygrid, kernel, np.sqrt(T * P * m), m, P, T
            )
            hy = np.cumsum((hy - hy[::-1]) * 0.5) * dqystep
            hz = histo_der(
                dq[:, 2::3], df[:, 2::3], dqzgrid, kernel, np.sqrt(T * P * m), m, P, T
            )
            hz = np.cumsum((hz - hz[::-1]) * 0.5) * dqzstep
            print(hx.sum() * dqxstep, hy.sum() * dqystep, hz.sum() * dqzstep)
            hxlist.append(hx)
            hylist.append(hy)
            hzlist.append(hz)

            # Computes the Fourier transform of the end to end vector.
            npx = hx * 0.0
            npy = hy * 0.0
            npz = hz * 0.0
            print("# computing FT for block #", x + 1)
            for i in range(len(hx)):
                npx[i] = (hx * np.cos(pxgrid[i] * dqxgrid)).sum() * dqxstep
                npy[i] = (hy * np.cos(pygrid[i] * dqygrid)).sum() * dqystep
                npz[i] = (hz * np.cos(pzgrid[i] * dqzgrid)).sum() * dqzstep

            nplistx.append(npx)
            nplisty.append(npy)
            nplistz.append(npz)

    # save the convoluted histograms of the end-to-end distances
    avghx = np.sum(np.asarray(hxlist), axis=0)
    errhx = np.std(np.asarray(hxlist), axis=0) / np.sqrt(n_block)
    avghy = np.sum(np.asarray(hylist), axis=0)
    errhy = np.std(np.asarray(hylist), axis=0) / np.sqrt(n_block)
    avghz = np.sum(np.asarray(hzlist), axis=0)
    errhz = np.std(np.asarray(hzlist), axis=0) / np.sqrt(n_block)

    norm_npx = avghx[(ns - 1) / 2]
    norm_npy = avghy[(ns - 1) / 2]
    norm_npz = avghz[(ns - 1) / 2]

    print("# Dx^2", np.dot(dqxgrid**2, avghx) * dqxstep / (bsize * n_block))
    print("# Dy^2", np.dot(dqygrid**2, avghy) * dqystep / (bsize * n_block))
    print("# Dz^2", np.dot(dqzgrid**2, avghz) * dqzstep / (bsize * n_block))

    print(
        "# px^2 (from the 2nd derivative of the histogram)",
        (
            30.0 * avghx[(ns - 1) / 2]
            - 16.0 * avghx[(ns - 1) / 2 + 1]
            - 16.0 * avghx[(ns - 1) / 2 - 1]
            + avghx[(ns - 1) / 2 - 2]
            + avghx[(ns - 1) / 2 + 2]
        )
        / dqxstep**2
        / norm_npx
        / 12.0,
    )
    print(
        "# py^2 (from the 2nd derivative of the histogram)",
        (
            30.0 * avghy[(ns - 1) / 2]
            - 16.0 * avghy[(ns - 1) / 2 + 1]
            - 16.0 * avghy[(ns - 1) / 2 - 1]
            + avghy[(ns - 1) / 2 - 2]
            + avghy[(ns - 1) / 2 + 2]
        )
        / dqystep**2
        / norm_npy
        / 12.0,
    )
    print(
        "# pz^2 (from the 2nd derivative of the histogram)",
        (
            30.0 * avghz[(ns - 1) / 2]
            - 16.0 * avghz[(ns - 1) / 2 + 1]
            - 16.0 * avghz[(ns - 1) / 2 - 1]
            + avghz[(ns - 1) / 2 - 2]
            + avghz[(ns - 1) / 2 + 2]
        )
        / dqzstep**2
        / norm_npz
        / 12.0,
    )

    np.savetxt("hxx.data", avghx)
    np.savetxt(
        str(prefix + "histo.data"),
        np.c_[
            dqxgrid,
            avghx / (bsize * n_block),
            errhx,
            dqygrid,
            avghy / (bsize * n_block),
            errhy,
            dqzgrid,
            avghz / (bsize * n_block),
            errhz,
        ],
    )

    # save the resulting momentum distribution for each axes
    avgnpx = np.sum(np.asarray(nplistx), axis=0) / 2.0 / np.pi / norm_npx
    avgnpy = np.sum(np.asarray(nplisty), axis=0) / 2.0 / np.pi / norm_npy
    avgnpz = np.sum(np.asarray(nplistz), axis=0) / 2.0 / np.pi / norm_npz

    errnpx = (
        np.std(np.asarray(nplistx), axis=0) * np.sqrt(n_block) / norm_npx / 2.0 / np.pi
    )
    errnpy = (
        np.std(np.asarray(nplisty), axis=0) * np.sqrt(n_block) / norm_npy / 2.0 / np.pi
    )
    errnpz = (
        np.std(np.asarray(nplistz), axis=0) * np.sqrt(n_block) / norm_npz / 2.0 / np.pi
    )

    np.savetxt(
        str(prefix + "np.data"),
        np.c_[pxgrid, avgnpx, errnpx, avgnpy, errnpy, avgnpz, errnpz],
    )

    # print the average value of p-square for each direction
    px2 = []
    py2 = []
    pz2 = []
    for i in range(n_block):
        px2.append(
            np.dot(pxgrid**2, np.asarray(nplistx)[i, :])
            * pxstep
            / norm_npx
            / 2.0
            / np.pi
            * n_block
        )
        py2.append(
            np.dot(pygrid**2, np.asarray(nplisty)[i, :])
            * pystep
            / norm_npy
            / 2.0
            / np.pi
            * n_block
        )
        pz2.append(
            np.dot(pzgrid**2, np.asarray(nplistz)[i, :])
            * pzstep
            / norm_npz
            / 2.0
            / np.pi
            * n_block
        )

    print("# number of blocks", n_block)
    print(
        "# px^2",
        np.mean(np.asarray(px2)),
        "sigmax",
        np.std(np.asarray(px2)) / np.sqrt(n_block),
    )
    print(
        "# py^2",
        np.mean(np.asarray(py2)),
        "sigmay",
        np.std(np.asarray(py2)) / np.sqrt(n_block),
    )
    print(
        "# pz^2",
        np.mean(np.asarray(pz2)),
        "sigmaz",
        np.std(np.asarray(pz2)) / np.sqrt(n_block),
    )
    print(
        "# p^2",
        np.mean(np.asarray(px2) + np.asarray(py2) + np.asarray(pz2)),
        "sigma",
        np.std(np.asarray(px2) + np.asarray(py2) + np.asarray(pz2)) / np.sqrt(n_block),
    )

    print("# time taken (s)", time.clock() - start)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument(
        "-qfile", type=str, default="", help="name of the bead positions file"
    )
    parser.add_argument(
        "-ffile", type=str, default="", help="name of the bead forces file"
    )
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
        "-ns",
        type=float,
        default=5000,
        help="Specify the number of point to use for the histogram",
    )
    parser.add_argument(
        "-skip", type=int, default=0, help="Specify the number of points to be skipped"
    )
    parser.add_argument(
        "-dint",
        type=float,
        default=10,
        help="Specify the positive extrema of the interval to build the histogram ([-dint,dint])",
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
        args.der,
        args.skip,
    )
