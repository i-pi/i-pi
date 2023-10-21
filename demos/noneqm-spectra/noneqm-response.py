#!/usr/bin/env python3
import numpy as np
import argparse
import xml.etree.ElementTree as et

def read_xml(file_in):
  tree = et.parse(file_in)
  root = tree.getroot()
  return tree, root

class Initializer():

    @staticmethod
    def load_from_xml(file_in, file_in_spec):
        """Loads i-pi input file `file_in` and spectra-specific input file `file_in_spec` """
        #----------------------------
        #Spectra-specific parameters.
        #----------------------------
        tree, root = read_xml(file_in_spec)
        #Number of steps in the correlation function/nonequilibrium dynamics.
        nsteps2 = int(root.find('./corr_steps').text)
        #Read operators.
        op = root.find('./op').text
        #Epsilon.
        epsilon = float(root.find('./epsilon').text)
        #Field polarization.
        field_pol = root.find('./field_pol').text
        #Output file name.
        try:
            out_name = root.find('./output').text + '_'
        except:
            out_name = ''

        #----------------------------
        #I-PI input parameters.
        #----------------------------
        tree, root = read_xml(file_in)
        #Predix for the i-pi output files.
        try:
            sim_name = root.find('./output').attrib['prefix'] # Read simulation name.
        except:
            sim_name = 'simulation'
        #Number of beads.
        try:
            nbeads = int(root.find('./system/initialize').attrib['nbeads']) # Read number of beads.
        except:
            nbeads = 1
        #Number of steps between nonequilibrium trajectories (taken from stride of chk files).
        step_traj = int(root.find('./output/checkpoint').attrib['stride'])
        #Stride for printing out dipole/polarizability values.
        for out in root.iter('trajectory'):
            if out.attrib['filename'] in op:
                step_print = int(out.attrib['stride'])
        #Total number of equilibrium steps.
        nsteps1 = int(root.find('./total_steps').text)
        #Size of dynamics timestep.
        delta_t = float(root.find('./system/motion/dynamics/timestep').text)
        #Ensemble temperature in K converted to beta in atomic units.
        beta = float(root.find('./system/ensemble/temperature').text)
        beta = 1.0 / (3.167e-6 * beta) #beta in atomic units.

        return Initializer(sim_name, nbeads, nsteps1//step_print, nsteps2//step_print, step_traj//step_print, step_print, op, epsilon, field_pol, beta, delta_t, out_name)

    def __init__(self, sim_name, nbeads, nsteps1, nsteps2, step1, step2, op, epsilon, field_pol, beta, delta_t, out_name):
        self.sim_name = sim_name
        self.nbeads = nbeads
        self.nsteps1 =  nsteps1
        self.nsteps2 = nsteps2
        self.step1 = step1
        self.step2 = step2
        self.op = op
        self.epsilon = epsilon
        self.field_pol = field_pol
        self.beta = beta
        self.delta_t = delta_t
        self.out_name = out_name

class ResponseFunction():

    """Implements two-dimensional equilibrium-nonequilibrium response function.
       The result is (beta / epsilon) <(op2_{+}(t2) - op2_{-}(t2)) op1(-t1)>, where op1 = op[0], op2 = op[1], 
       and op is given in input and has two components that are either dip (dipole) or pol (poalrizability).
       Dipole and polarizability are projected on the field polarization given by the field_pol parameter.
       Output is related to the response function of two-dimensional IR-Raman spectroscopy.
    """

    def __init__(self, input_param):
        self.sim_name = input_param.sim_name
        self.out_name = input_param.out_name
        self.nbeads = input_param.nbeads
        self.op = input_param.op
        self.nsteps1 = input_param.nsteps1
        self.nsteps2 = input_param.nsteps2
        self.step1 = input_param.step1
        self.step2 = input_param.step2

        self.fmt_bead = (
            "{0:0"
            + str(int(1 + np.floor(np.log(self.nbeads) / np.log(10))))
            + "d}"
        )
        self.epsilon = input_param.epsilon
        self.beta = input_param.beta
        self.delta_t = input_param.delta_t
        self.field_pol = [ int(i) for i in input_param.field_pol.split(',') ]
        self.op = [ i.strip() for i in self.op.split(',') ]

    def process(self):
        #Get values of dip/pol along field_pol direction and average over beads. 
        x = []
        op_index = 0
        for op, field_pol in zip(self.op, self.field_pol):
            x_single = 0
            for bead in range(self.nbeads):
                if op_index == 0:
                    x_single += np.loadtxt(self.sim_name + '.' + op + '_' + self.fmt_bead.format(bead))[:, field_pol]
                else:
                    x_single += np.loadtxt(self.sim_name + '_noneqm.' + op + '_' + self.fmt_bead.format(bead))[:, field_pol]
            op_index += 1
            x_single /= self.nbeads
            x.append(x_single)

        #Compute result:
        c = 0 # Count nonequilibrium events.
        corr = 0
        i_step = self.nsteps2 + 1 #Number of steps for t1 and t2.
        for i in range(self.step1 + 1, len(x[0])+1, self.step1): #Start at the first available checkpoint from eq dynamics.
            if i >= i_step: #Check if we have enough dynamics to go backward -t1.
                #Backward equilibrium trajectory, first operator.
                eq = np.flip(x[0][i - i_step : i]) #Equilibrium dynamics part from time t (index i) to t-t1.
                #Difference between forward nonequilibrium trajectories, second operator.
                neq = x[1][2 * c * i_step : (2 * c + 1) * i_step] - x[1][(2 * c + 1) * i_step : 2 * (c + 1) * i_step] #Neq dynamics.
                #Compute corr[j, k] = eq[j] * neq[k].
                corr+=np.outer(eq, neq)
                c+=1
        corr *= self.beta / (c * self.epsilon)
        corr = np.gradient(corr, self.step2 * self.delta_t, axis=0, edge_order=2) 

        np.savetxt(self.out_name + self.op[0] + '_' + self.op[1] + '_neq_2d.dat', corr)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Equilibrium-nonequilibrium response functions.')
    parser.add_argument('--input_ipi', default='input.xml', help='I-PI input xml file.')
    parser.add_argument('--input_spec', default='input_spec.xml', help='Response function input xml file.')
    args = parser.parse_args()

    init_param = Initializer.load_from_xml(args.input_ipi, args.input_spec)
    process = ResponseFunction(init_param)
    process.process()
