"""Contains all methods to evalaute potential energy and forces for indistinguishable particles.
Used in /engine/normalmodes.py
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import math
import numpy as np
from ipi.utils.depend import *
import ipi.utils.units as units
import sys


def kth_diag_indices(a, k):
    """
    Indices to access matrix k-diagonals in numpy.
    https://stackoverflow.com/questions/10925671/numpy-k-th-diagonal-indices
    """
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols


class ExchangePotential:
    def __init__(self, nbeads, nbosons, bead_mass):
        assert nbosons > 0
        self._N = nbosons
        self._P = nbeads
        self._particle_mass = bead_mass

        self._bead_diff_intra = None
        self._bead_diff_inter_first_last_bead = None
        self._E_from_to = None
        self._V = None
        self._V_backward = None

        # creates a placeholder for omegan2
        self.omegan2 = depend_value("omegan2", value=0.0)

    
    def bind(self, beads, ensemble):
        self.beads = beads
        self.nm = normalmodes
        self.ensemble = ensemble
        self.nbeads -

        self.betaP = depend_value("betaP",
                                   func=lambda: 1.0 / (self.beads.nbeads * units.Constants.kb * self.ensemble.temp),
                                   dependencies = [self.ensemble._temp]
                                   )
        
        dpipe(self.nm._omegan2, self._omegan2)

        self._V = depend_value("V", func=self._evaluate_prefix_V,
                               dependencies=[self.beads._q, self._betaP])
        
        self._vspring_and_fspring = depend_value("vspring_and_fspring", 
                                                 func=self.get_vspring_and_fspring,
                                                 dependencies=[]
                                                 )

    def set_coordinates(self, q):
        self._q = q

        # self._bead_diff_intra[j] = [r^{j+1}_0 - r^{j}_0, ..., r^{j+1}_{N-1} - r^{j}_{N-1}]
        self._bead_diff_intra = np.diff(self._q, axis=0)
        # self._bead_dist_inter_first_last_bead[l][m] = r^0_{l} - r^{P-1}_{m}
        self._bead_diff_inter_first_last_bead = (
            self._q[0, :, np.newaxis, :] - self._q[self._P - 1, np.newaxis, :, :]
        )

        # cycle energies:
        # self._E_from_to[u, v] is the ring polymer energy of the cycle on particle indices u,...,v
        self._E_from_to = self._evaluate_cycle_energies()

        # prefix potentials:
        # self._V[l] = V^[1,l+1]
        self._V = self._evaluate_prefix_V()

        # suffix potentials:
        # self._V_backward[l] = V^[l+1, N]
        self._V_backward = self._evaluate_suffix_V()

    def get_vspring_and_fspring(self):
        """
        Returns spring potential and forces for bosons.
        """
        F = self.evaluate_spring_forces()

        return [self.V_all(), F]

    def V_all(self):
        """
        Returns the potential on all particles: V^[1,N]
        """
        return self._V[self._N]

    def _spring_potential_prefix(self):
        """
        Helper function: the term for the spring constant as used in spring potential expressions
        """
        return 0.5 * self._particle_mass * self.omegan2

    def _spring_force_prefix(self):
        """
        Helper function: the term for the spring constant as used in spring force expressions
        """
        return (-1.0) * self._particle_mass * self.omegan2

    def _evaluate_cycle_energies(self):
        """
        Evaluate all the cycle energies, as outlined in Eqs. 5-7 of arXiv:2305.18025.
        Returns an upper-triangular matrix, Emks[u,v] is the ring polymer energy of the cycle on u,...,v.
        """
        # using column-major (Fortran order) because uses of the array are mostly within the same column
        Emks = np.zeros((self._N, self._N), dtype=float, order="F")

        intra_spring_energies = np.sum(self._bead_diff_intra**2, axis=(0, -1))
        spring_energy_first_last_bead_array = np.sum(
            self._bead_diff_inter_first_last_bead**2, axis=-1
        )

        # for m in range(self._N):
        #     Emks[m][m] = self._spring_potential_prefix() * \
        #                     (intra_spring_energies[m] + spring_energy_first_last_bead_array[m, m])
        Emks[np.diag_indices_from(Emks)] = self._spring_potential_prefix() * (
            intra_spring_energies + np.diagonal(spring_energy_first_last_bead_array)
        )

        for s in range(self._N - 1 - 1, -1, -1):
            # for m in range(s + 1, self._N):
            #     Emks[s][m] = Emks[s + 1][m] + self._spring_potential_prefix() * (
            #             - spring_energy_first_last_bead_array[s + 1, m]
            #             + intra_spring_energies[s]
            #             + spring_energy_first_last_bead_array[s + 1, s]
            #             + spring_energy_first_last_bead_array[s, m])
            Emks[s, (s + 1) :] = Emks[
                s + 1, (s + 1) :
            ] + self._spring_potential_prefix() * (
                -spring_energy_first_last_bead_array[s + 1, (s + 1) :]
                + intra_spring_energies[s]
                + spring_energy_first_last_bead_array[s + 1, s]
                + spring_energy_first_last_bead_array[s, (s + 1) :]
            )

        return Emks

    def _evaluate_prefix_V(self):
        """
        Evaluate V^[1,1], V^[1,2], ..., V^[1,N], as outlined in Eq. 3 of arXiv:2305.18025.
        (In the code, particle indices start from 0.)
        Returns a vector of these potentials, in this order.
        Assumes that the cycle energies self._E_from_to have been computed.
        """
        V = np.zeros(self._N + 1, float)

        for m in range(1, self._N + 1):
            # For numerical stability
            subdivision_potentials = V[:m] + self._E_from_to[:m, m - 1]
            Elong = np.min(subdivision_potentials)

            # sig = 0.0
            # for u in range(m):
            #   sig += np.exp(- self.betaP *
            #                (V[u] + self._E_from_to[u, m - 1] - Elong) # V until u-1, then cycle from u to m
            #                 )
            sig = np.sum(np.exp(-self.betaP * (subdivision_potentials - Elong)))
            assert sig != 0.0 and np.isfinite(sig)
            V[m] = Elong - np.log(sig / m) / self.betaP

        return V

    def _evaluate_suffix_V(self):
        """
        Evaluate V^[1,N], V^[2,N], ..., V^[N,N], as outlined in Eq. 16 of arXiv:2305.18025.
        (In the code, particle indices start from 0.)
        Returns a vector of these potentials, in this order.
        Assumes that both the cycle energies self._E_from_to and prefix potentials self._V have been computed.
        """
        RV = np.zeros(self._N + 1, float)

        for l in range(self._N - 1, 0, -1):
            # For numerical stability
            subdivision_potentials = self._E_from_to[l, l:] + RV[l + 1 :]
            Elong = np.min(subdivision_potentials)

            # sig = 0.0
            # for p in range(l, self._N):
            #     sig += 1 / (p + 1) * np.exp(- self.betaP * (self._E_from_to[l, p] + RV[p + 1]
            #                                                 - Elong))
            sig = np.sum(
                np.reciprocal(np.arange(l + 1.0, self._N + 1))
                * np.exp(-self.betaP * (subdivision_potentials - Elong))
            )
            assert sig != 0.0 and np.isfinite(sig)
            RV[l] = Elong - np.log(sig) / self.betaP

        # V^[1,N]
        RV[0] = self._V[-1]

        return RV

    def evaluate_spring_forces(self):
        """
        Evaluate the ring polymer forces on all the beads, as outlined in Eq. 13, 17-18 of arXiv:2305.18025.
        (In the code, particle indices start from 0.)
        Returns an array, F[j, l, :] is the force (3d vector) on bead j of particle l.
        Assumes that both the cycle energies self._E_from_to, the prefix potentials self._V,
        and the suffix potentials self._V_backward have been computed.
        """
        F = np.zeros((self._P, self._N, 3), float)

        # force on intermediate beads
        #
        # for j in range(1, self._P - 1):
        #   for l in range(self._N):
        #         F[j, l, :] = self._spring_force_prefix() * (-self._bead_diff_intra[j][l] +
        #                                                     self._bead_diff_intra[j - 1][l])
        F[1:-1, :, :] = self._spring_force_prefix() * (
            -self._bead_diff_intra[1:, :]
            + np.roll(self._bead_diff_intra, axis=0, shift=1)[1:, :]
        )

        # force on endpoint beads
        #
        connection_probs = np.zeros((self._N, self._N), float)
        # close cycle probabilities:
        # for u in range(0, self._N):
        #     for l in range(u, self._N):
        #         connection_probs[l][u] = 1 / (l + 1) * \
        #                np.exp(- self.betaP *
        #                        (self._V[u] + self._E_from_to[u, l] + self._V_backward[l+1]
        #                         - self.V_all()))
        tril_indices = np.tril_indices(self._N, k=0)
        connection_probs[tril_indices] = (
            # np.asarray([1 / (l + 1) for l in range(self._N)])[:, np.newaxis] *
            np.reciprocal(np.arange(1.0, self._N + 1))[:, np.newaxis]
            * np.exp(
                -self.betaP
                * (
                    # np.asarray([self._V(u - 1) for u in range(self._N)])[np.newaxis, :]
                    self._V[np.newaxis, :-1]
                    # + np.asarray([(self._E_from_to[u, l] if l >= u else 0) for l in range(self._N)
                    #                   for u in range(self._N)]).reshape((self._N, self._N))
                    + self._E_from_to.T
                    # + np.asarray([self._V_backward(l + 1) for l in range(self._N)])[:, np.newaxis]
                    + self._V_backward[1:, np.newaxis]
                    - self.V_all()
                )
            )
        )[tril_indices]

        # direct link probabilities:
        # for l in range(self._N - 1):
        #     connection_probs[l][l+1] = 1 - (np.exp(- self.betaP * (self._V[l + 1] + self._V_backward[l + 1] -
        #                                         self.V_all())))
        superdiagonal_indices = kth_diag_indices(connection_probs, k=1)
        connection_probs[superdiagonal_indices] = 1 - (
            np.exp(
                -self.betaP * (self._V[1:-1] + self._V_backward[1:-1] - self.V_all())
            )
        )

        # on the last bead:
        #
        # for l in range(self._N):
        #     force_from_neighbor = np.empty((self._N, 3))
        #     for next_l in range(max(l + 1, self._N)):
        #         force_from_neighbor[next_l, :] = self._spring_force_prefix() * \
        #                         (-self._bead_diff_inter_first_last_bead[next_l, l] + self._bead_diff_intra[-1, l])
        #     F[-1, l, :] = sum(connection_probs[l][next_l] * force_from_neighbor[next_l]
        #                       for next_l in range(self._N))
        #
        # First vectorization:
        # for l in range(self._N):
        #     force_from_neighbors = np.empty((self._N, 3))
        #     force_from_neighbors[:, :] = self._spring_force_prefix() * \
        #                         (-self._bead_diff_inter_first_last_bead[:, l] + self._bead_diff_intra[-1, l])
        #     F[-1, l, :] = np.dot(connection_probs[l][:], force_from_neighbors)
        force_from_neighbors = self._spring_force_prefix() * (
            -np.transpose(self._bead_diff_inter_first_last_bead, axes=(1, 0, 2))
            + self._bead_diff_intra[-1, :, np.newaxis]
        )
        # F[-1, l, k] = sum_{j}{force_from_neighbors[l][j][k] * connection_probs[l,j]}
        F[-1, :, :] = np.einsum("ljk,lj->lk", force_from_neighbors, connection_probs)

        # on the first bead:
        #
        # for l in range(self._N):
        #     force_from_neighbor = np.empty((self._N, 3))
        #     for prev_l in range(l - 1, self._N):
        #         force_from_neighbor[prev_l, :] = self._spring_force_prefix() * \
        #                            (-self._bead_diff_intra[0, l] + self._bead_diff_inter_first_last_bead[l, prev_l])
        #     F[0, l, :] = sum(connection_probs[prev_l][l] * force_from_neighbor[prev_l]
        #                      for prev_l in range(self._N))
        #
        # First vectorization:
        #
        # for l in range(self._N):
        #     force_from_neighbors = np.empty((self._N, 3))
        #     force_from_neighbors[:, :] = self._spring_force_prefix() * \
        #                              (-self._bead_diff_intra[0, l] + self._bead_diff_inter_first_last_bead[l, :])
        #     F[0, l, :] = np.dot(connection_probs[:, l], force_from_neighbors)
        #
        force_from_neighbors = self._spring_force_prefix() * (
            self._bead_diff_inter_first_last_bead
            - self._bead_diff_intra[0, :, np.newaxis]
        )
        # F[0, l, k] = sum_{j}{force_from_neighbors[l][j][k] * connection_probs[j,l]}
        F[0, :, :] = np.einsum("ljk,jl->lk", force_from_neighbors, connection_probs)

        return F

    def get_distinct_probability(self):
        """
        Evaluate the probability of the configuration where all the particles are separate.
        """
        return np.exp(
            -self.betaP * (np.trace(self._E_from_to) - self.V_all())
            - math.log(np.math.factorial(self._N))  # (1.0 / np.math.factorial(self._N))
        )

    def get_longest_probability(self):
        """
        Evaluate the probability of a configuration where all the particles are connected,
        divided by 1/N. Notice that there are (N-1)! permutations of this topology
        (all represented by the cycle 0,1,...,N-1,0); this cancels the division by 1/N.
        """
        return np.exp(-self.betaP * (self._E_from_to[0, -1] - self.V_all()))

    def get_kinetic_td(self):
        """Implementation of the Hirshberg-Rizzi-Parrinello primitive
        kinetic energy estimator for identical particles.
        Corresponds to Eqns. (4)-(5) in SI of pnas.1913365116.
        """
        est = np.zeros(self._N + 1)

        for m in range(1, self._N + 1):
            sig = 0.0

            # Numerical stability - Xiong-Xiong method (arXiv.2206.08341)
            e_tilde = sys.float_info.max
            for k in range(m, 0, -1):
                e_tilde = min(e_tilde, self._E_from_to[m - k, m - 1] + self._V[m - k])

            for k in range(m, 0, -1):
                E_kn_val = self._E_from_to[m - k, m - 1]
                sig += (est[m - k] - E_kn_val) * np.exp(
                    -self.betaP * (E_kn_val + self._V[m - k] - e_tilde)
                )

            sig_denom_m = m * np.exp(-self.betaP * (self._V[m] - e_tilde))

            est[m] = sig / sig_denom_m

        factor = 1.5 * self._N / self.betaP

        return factor + est[self._N] / self._P

    def get_fermionic_sign(self):
        """
        The average permutation sign as defined in Eq. (9) https://doi.org/10.1063/5.0008720,
        which can be used to reweight observables to obtain fermionic statistics.
        """
        return self._get_fermionic_potential_exp() / np.exp(-self.betaP * self._V[-1])

    def _get_fermionic_potential_exp(self):
        """
        Exponential of the fermionic pseudo-potential defined by
        the recurrence relation in Eq. (5) of https://doi.org/10.1063/5.0008720.
        Numerically unstable since it does not use log-sum-exp trick, seeing that the
        sum of exponentials could be negative.
        """
        xi = -1
        W = np.empty(self._N + 1, float)
        W[0] = 1.0

        for m in range(1, self._N + 1):
            perm_sign = np.array([xi ** (k - 1) for k in range(m, 0, -1)])
            W[m] = (1.0 / m) * np.sum(
                perm_sign * W[:m] * np.exp(-self.betaP * self._E_from_to[:m, m - 1])
            )

        return W[-1]
