"""
reweight object
"""
import sys

import numpy as np
import time


def fancy_print(string):
    """
    create a fancy print
    :param string: the string for printing out
    :return:
    """
    print("#Reweight: {0}".format(string))

class Reweight():
    au2eV = 2.72113838565563E+01
    def __init__(self, config_file):
        # load json file
        import json
        with open(config_file, 'r') as f:
            config = json.load(f)
        # load the pot file
        self._frame_pot_list = np.loadtxt(config["pot_file"])
        fancy_print("Load potential data file {0}".format(config["pot_file"]))
        self._quant_list = np.loadtxt(config["quant_file"])
        fancy_print("Load quantity data file {0}".format(config["quant_file"]))
        #special case for 1 d array with other shape
        if self._quant_list.ndim == 1 :
            self._quant_list = np.expand_dims(self._quant_list, axis=0)
        else:
            self._quant_list = self._quant_list.T

        self._num_quant = self._quant_list.shape[0]
        quantity_frame = self._quant_list.shape[1]
        fancy_print("Read Quantity Number: {0}".format(self._num_quant))
        # reading the neccessary parameter
        self._method = config["reweight_method"]
        fancy_print("Reweight Method: {0}".format(self._method))
        self._temperature = float(config["temp"])
        fancy_print("Read Temperature: {0} K".format(self._temperature))
        self._num_frame = self._frame_pot_list.shape[0]
        fancy_print("Read Frame Number: {0}".format(self._num_frame))
        self._num_pot = self._frame_pot_list.shape[1]
        fancy_print("Read Machine Learning Potential Number: {0}".format(self._num_pot))
        # reading the old and new alpha value
        self._old_alpha = config["old_alpha"]
        fancy_print("Read old alpha value: {0}".format(self._old_alpha))
        self._new_alpha = config["new_alpha"]
        fancy_print("Read new alpha value: {0}".format(self._new_alpha))
        # check the frame equality in pot file and quant file
        if self._num_frame != quantity_frame:
            fancy_print("the number of frame in potential file "
                        "is not equal to "
                        "the number of frame in quantity file")
            sys.exit("exit the program")

        # if additional weight parameter is added
        if "extra_weight" in config:
            fancy_print("a extra weight file is detected")
            self._extra_weight = np.loadtxt(config["extra_weight"])
            self._has_extra_file = True
            fancy_print("Load extra weight file {0}".format(config["extra_weight"]))
        else:
            self._has_extra_file = False



    def run(self):
        """
        main program run the reweight
        :return:
        """
        self.compute_ave_pot()
        self.compute_rescaled_pot()
        if self._method == "Direct":
            self.compute_weight_matrix()
            self.direct_compute_reweight_quantity(self._has_extra_file)
        elif self._method == "Asymptotic":
            self.asymptotic_compute_reweight_quantity(self._has_extra_file)
        else:
            fancy_print("No Implemented Method!!!")
            sys.exit("exiting")

    def compute_ave_pot(self):
        """
        compute the average potential
        :return:
        """
        fancy_print("Starting compute the average potential")
        start = time.time()
        self._frame_ave_pot = self._frame_pot_list.mean(axis=1)
        end = time.time()
        fancy_print("Computing average potential finished, "
                    "time: {0:.6f}s".format(end-start))

    def compute_rescaled_pot(self):
        """
        compute the rescaled potential with from old alpha to new alpha
        """
        fancy_print("Starting compute the rescaled potential")
        start = time.time()

        self._delta_pot = np.zeros([self._num_frame, self._num_pot])
        for i in range(self._num_frame):
            for j in range(self._num_pot):
                self._delta_pot[i][j] = self._frame_pot_list[i][j] - self._frame_ave_pot[i]

        self._delta_pot = self._new_alpha / self._old_alpha * self._delta_pot
        for i in range(self._num_frame):
            for j in range(self._num_pot):
                self._frame_pot_list[i][j] = self._delta_pot[i][j] + self._frame_ave_pot[i]

        end = time.time()
        fancy_print("Computing rescaled potential finished, "
                    "time: {0:.6f}s".format(end-start))



    def compute_weight_matrix(self):
        kB = 8.617333262145e-5
        beta = 1.0 / (kB * self._temperature)
        self._weight_matrix = np.zeros([self._num_frame, self._num_pot])
        fancy_print("Starting compute the weight matrix")
        start = time.time()
        for i in range(self._num_frame):
            for j in range(self._num_pot):
                self._weight_matrix[i][j] = (self._frame_pot_list[i][j] - self._frame_ave_pot[i])
        self._weight_matrix = np.exp(self._weight_matrix * (-1.0 * beta))
        #normalize the weight_matrix
        norm = self._weight_matrix.sum(axis=0)
        for i in range(self._num_pot):
            self._weight_matrix[:, i] /= norm[i]
        end = time.time()
        fancy_print("Computing weight matrix finished, "
                    "time: {0:.6f}s".format(end-start))
        np.savetxt("weight_matrix.dat", self._weight_matrix)
        fancy_print("Save weight to file: weight_matrix.dat")

    def direct_compute_reweight_quantity(self, has_extra_file):
        fancy_print("Starting compute the reweight quantity")
        start = time.time()
        if has_extra_file:
            self._weight_matrix *= self._extra_weight
        reweight_quant = np.matmul(self._quant_list, self._weight_matrix)
        end = time.time()
        fancy_print("Computing reweight quantity finished, "
                    "time: {0:.6f}s".format(end-start))

        np.savetxt("direct_reweight_quantity.dat", reweight_quant)
        fancy_print("Save reweight quantity to file: "
                    "direct_reweight_quantity.dat")

    def asymptotic_compute_reweight_quantity(self, has_extra_file):
        """
        reweight the quantity by asymptotic approximation
        ref:
        Michele Ceriotti,  Guy A. R. Brain, Oliver Riordan and David E.
        Manolopoulos 2011The inefficiency of re-weighted sampling and
        the curse of system size in high-order path integration
        Proc. R. Soc. A.4682–17
        http://doi.org/10.1098/rspa.2011.0413
        :return:
        """
        kB = 8.617333262145e-5
        beta = 1.0 / (kB * self._temperature)

        fancy_print("Starting compute the reweight quantity")
        start = time.time()

        # get the h matrix h = beta(E_i - E_ave)
        self._h_matrix = np.zeros([self._num_frame, self._num_pot])
        for i in range(self._num_frame):
            for j in range(self._num_pot):
                self._h_matrix[i][j] = (self._frame_pot_list[i][j] - self._frame_ave_pot[i])*beta
        if has_extra_file:
            pass
            #self._h_matrix -= np.log(self._extra_weight)

        # get h_ave
        self._h_ave = self._h_matrix.mean(axis=0)
        # get quantity average, i.e. <q>
        self._quant_ave = self._quant_list.mean(axis=1)
        #get <qh>
        self._qh_ave = np.matmul(self._quant_list, self._h_matrix)/self._num_frame
        #get the asymptotic quantity
        asym_reweight_quant = np.zeros([self._num_quant, self._num_pot])
        for i in range(self._num_quant):
            for j in range(self._num_pot):
                asym_reweight_quant[i][j] = self._quant_ave[i] - self._qh_ave[i][j] \
                                            + self._quant_ave[i]*self._h_ave[j]

        end = time.time()
        fancy_print("Computing reweight quantity finished, "
                    "time: {0:.6f}s".format(end - start))
        np.savetxt("asymptotic_reweight_quantity.dat", asym_reweight_quant)
        fancy_print("Save reweight quantity to file: "
                    "asymptotic_reweight_quantity.dat")











