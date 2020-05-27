"""Minimal compressed sparse row and compressed sparse column sparse matrix implementation,
used only for storage and retrieval but (ATM) not for calculations.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2019 i-PI developers
# See the "licenses" directory for full license information.

import numpy as np

##########################################################################################

# Load sparse matrix file in binary format
# Pass a function capable of loading the file to avoid np.load (loader='numpy')


def load_npz(file, loader="numpy"):
    if loader == "numpy":
        data = np.load(file)
    else:
        data = loader(file)
    if data["kind"] == "csr":
        return csr_matrix(
            a=data["a"], ia=data["ia"], ja=data["ja"], m=data["m"], n=data["n"]
        )
    else:
        return csc_matrix(
            a=data["a"], ia=data["ia"], ja=data["ja"], m=data["m"], n=data["n"]
        )


##########################################################################################

# Save sparse matrix file in binary format
# Pass a function capable of saving the file to avoid np.savez (saver='numpy')


def save_npz(file, matrix, saver="numpy", compressed=True):
    a = matrix.a
    ia = matrix.ia
    ja = matrix.ja
    m = matrix.m
    n = matrix.n
    kind = matrix.kind
    if saver == "numpy":
        if compressed is True:
            np.savez_compressed(file, a=a, ia=ia, ja=ja, m=m, n=n, kind=kind)
        else:
            np.savez(file, a=a, ia=ia, ja=ja, m=m, n=n, kind=kind)
    else:
        saver(file, compressed=compressed, a=a, ia=ia, ja=ja, m=m, n=n, kind=kind)


##########################################################################################

# Generic base class


class sparse_matrix(object):
    def __add__(self, b):
        densea = self.toarray()
        denseb = b.toarray()
        return csr_matrix(densea + denseb)

    def __sub__(self, b):
        densea = self.toarray()
        denseb = b.toarray()
        return csr_matrix(densea - denseb)

    def __mul__(self, b):
        densea = self.toarray()
        denseb = b.toarray()
        return csr_matrix(densea * denseb)

    # dot product
    def dot(self, b):
        densea = self.toarray()
        denseb = b.toarray()
        return csr_matrix(np.dot(densea, denseb))

    # fraction of non-zero elements in matrix
    def density(self):
        return float(len(self.a)) / float(self.m * self.n)


##########################################################################################

# Compressed sparse row format


class csr_matrix(sparse_matrix):
    def __init__(self, nparray=None, **kwargs):
        self.kind = "csr"
        # Construct by converting a numpy array to csr
        if "a" not in kwargs:
            self.m = nparray.shape[0]
            self.n = nparray.shape[1]
            nonzeroids = np.where(nparray != 0)
            self.a = nparray[nonzeroids]
            self.ia = np.zeros(self.m + 1, dtype=np.int64)
            for row in nonzeroids[0]:
                self.ia[row + 1] += 1
            self.ia = np.cumsum(self.ia)
            self.ja = nonzeroids[1]
        # Construct from sparse format arrays
        else:
            self.a = kwargs["a"]
            self.ia = kwargs["ia"]
            self.ja = kwargs["ja"]
            self.m = kwargs["m"]
            self.n = kwargs["n"]

    # Convert to dense format
    def toarray(self):
        nparray = np.zeros((self.m, self.n), dtype=self.a.dtype)
        for row in range(self.m):
            nparray[row, self.ja[self.ia[row] : self.ia[row + 1]]] = self.a[
                self.ia[row] : self.ia[row + 1]
            ]
        return nparray


# #########################################################################################

# Compressed sparse column format


class csc_matrix(sparse_matrix):
    def __init__(self, nparray=None, **kwargs):
        self.kind = "csc"
        # Construct by converting a numpy array to csr
        if "a" not in kwargs:
            self.m = nparray.shape[0]
            self.n = nparray.shape[1]
            nonzeroids = np.where(nparray.T != 0)
            self.a = nparray.T[nonzeroids]
            self.ia = np.zeros(self.n + 1, dtype=np.int64)
            for col in nonzeroids[0]:
                self.ia[col + 1] += 1
            self.ia = np.cumsum(self.ia)
            self.ja = nonzeroids[1]
        # Construct from sparse format arrays
        else:
            self.a = kwargs["a"]
            self.ia = kwargs["ia"]
            self.ja = kwargs["ja"]
            self.m = kwargs["m"]
            self.n = kwargs["n"]

    # Convert to dense format
    def toarray(self):
        nparray = np.zeros((self.m, self.n), dtype=self.a.dtype)
        for col in range(self.n):
            nparray[self.ja[self.ia[col] : self.ia[col + 1]], col] = self.a[
                self.ia[col] : self.ia[col + 1]
            ]
        return nparray


# #########################################################################################
